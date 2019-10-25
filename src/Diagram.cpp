#include "Diagram.hpp"

#include "KahanAccumulator.hpp"
#include "timings.hpp"

#include <omp.h>
#include <boost/range/adaptor/indexed.hpp>

Complex resolve_request1(std::vector<TraceRequest> const &trace_requests,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  TimingScope<1> timing_scope("resolve_request1");

  assert(ssize(trace_requests) == 1);
  auto const &trace_request0 = trace_requests.at(0);
  auto const &locations0 = trace_request0.locations;

  if (q.trace_factories.count(trace_request0.tr_name) == 0) {
    std::cerr << trace_request0.tr_name << std::endl;
  }
  auto const &x0 = q.trace_factories.at(trace_request0.tr_name)
                       ->get(slice_pair, locations0)
                       .at(trace_request0.tr_id);

  return std::accumulate(std::begin(x0), std::end(x0), Accumulator<Complex>{}).value() /
         static_cast<double>(x0.size());
}

Complex resolve_request2(std::vector<TraceRequest> const &trace_requests,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  TimingScope<1> timing_scope("resolve_request2");

  assert(ssize(trace_requests) == 2);
  auto const &trace_request0 = trace_requests.at(0);
  auto const &locations0 = trace_request0.locations;
  auto const &x0 = q.trace_factories.at(trace_request0.tr_name)
                       ->get(slice_pair, locations0)
                       .at(trace_request0.tr_id);
  DilutedTraces t0{x0, false};

  auto const &trace_request1 = trace_requests.at(1);
  auto const &locations1 = trace_request1.locations;
  auto const &x1 = q.trace_factories.at(trace_request1.tr_name)
                       ->get(slice_pair, locations1)
                       .at(trace_request1.tr_id);
  DilutedTraces t1{x1, false};

  return inner_product(t0, t1);
}

Complex resolve_request3(std::vector<TraceRequest> const &trace_requests,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  TimingScope<1> timing_scope("resolve_request3");

  assert(ssize(trace_requests) == 2);
  auto const &trace_request0 = trace_requests.at(0);
  auto const &locations0 = trace_request0.locations;
  auto const &x0 = q.trace_factories.at(trace_request0.tr_name)
                       ->get(slice_pair, locations0)
                       .at(trace_request0.tr_id);
  DilutedTraces t0{x0, false};

  auto const &trace_request1 = trace_requests.at(1);
  auto const &locations1 = trace_request1.locations;
  auto const &x1 = q.trace_factories.at(trace_request1.tr_name)
                       ->get(slice_pair, locations1)
                       .at(trace_request1.tr_id);
  DilutedTraces t1{x1, false};

  auto const &trace_request2 = trace_requests.at(2);
  auto const &locations2 = trace_request2.locations;
  auto const &x2 = q.trace_factories.at(trace_request2.tr_name)
                       ->get(slice_pair, locations2)
                       .at(trace_request2.tr_id);
  DilutedTraces t2{x2, false};

  return inner_product(t0, t1, t2);
}

Complex resolve_request(std::vector<TraceRequest> const &trace_requests,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  if (ssize(trace_requests) == 1) {
    return resolve_request1(trace_requests, slice_pair, q);
  } else if (ssize(trace_requests) == 2) {
    return resolve_request2(trace_requests, slice_pair, q);
  } else if (ssize(trace_requests) == 3) {
    return resolve_request3(trace_requests, slice_pair, q);
  } else {
    throw std::runtime_error("This many traces are not implemented yet.");
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
