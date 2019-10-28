#include "Diagram.hpp"

#include "KahanAccumulator.hpp"
#include "timings.hpp"

#include <omp.h>
#include <boost/range/adaptor/indexed.hpp>

void request_request(std::vector<TraceRequest> const &trace_requests,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) {
  TimingScope<1> timing_scope("request_request");

  for (auto const &trace_request : trace_requests) {
    q.trace_factories.at(trace_request.tr_name)
        ->request(slice_pair, trace_request.locations);
  }
}

Complex resolve_request(std::vector<TraceRequest> const &trace_requests,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  TimingScope<1> timing_scope("resolve_request");

  std::vector<DilutedTraces> dt;
  dt.reserve(trace_requests.size());

  for (auto const &trace_request : trace_requests) {
    auto const &locations = trace_request.locations;
    auto const &x = q.trace_factories.at(trace_request.tr_name)
                        ->get(slice_pair, locations)
                        .at(trace_request.tr_id);
    DilutedTraces t{x, false};
    dt.push_back(t);
  }

  if (ssize(trace_requests) == 1) {
    return std::accumulate(
               std::begin(dt[0].traces), std::end(dt[0].traces), Accumulator<Complex>{})
               .value() /
           static_cast<double>(dt[0].traces.size());
  } else if (ssize(trace_requests) == 2) {
    return inner_product(dt[0], dt[1]);
  } else if (ssize(trace_requests) == 3) {
    return inner_product(dt[0], dt[1], dt[2]);
  } else {
    throw std::runtime_error("This many traces are not implemented yet.");
  }
}

void Diagram::request(int const t, BlockIterator const &slice_pair, DiagramParts &q) {
  TimingScope<1> timing_scope("Diagram::request", name());

  request_impl(t, slice_pair, q);
}

void Diagram::request_impl(int const t,
                           BlockIterator const &slice_pair,
                           DiagramParts &q) {
  TimingScope<1> timing_scope("Diagram::request_impl", name());

  for (int i = 0; i != ssize(correlator_requests()); ++i) {
    auto const &request = correlator_requests()[i];
    request_request(request.trace_requests, slice_pair, q);
  }
}

void Diagram::assemble(int const t, BlockIterator const &slice_pair, DiagramParts &q) {
  TimingScope<1> timing_scope("Diagram::assemble", name());

  assemble_impl(t, slice_pair, q);
}

void Diagram::assemble_impl(int const t,
                            BlockIterator const &slice_pair,
                            DiagramParts &q) {
  TimingScope<1> timing_scope("Diagram::assemble_impl", name());

  for (int i = 0; i != ssize(correlator_requests()); ++i) {
    auto const &request = correlator_requests()[i];
    auto const &number = resolve_request(request.trace_requests, slice_pair, q);
    correlator_.at(t).at(i) += number;
  }
}

void Diagram::write() {
  TimingScope<1> timing_scope("Diagram::write", name());

  WriteHDF5Correlator filehandle(
      output_path_, name(), output_filename_, make_comp_type<Complex>());

  std::vector<Complex> one_corr(Lt_);

  for (int i = 0; i != ssize(correlator_requests()); ++i) {
    for (int t = 0; t < Lt_; ++t) {
      one_corr[t] = correlator_[t][i].value() / static_cast<double>(Lt_);
    }
    filehandle.write(one_corr, correlator_requests()[i].hdf5_dataset_name);
  }
}
