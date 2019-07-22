#pragma once

#include "Correlators.hpp"
#include "DilutedFactorFactory.hpp"
#include "DilutedProductFactory.hpp"
#include "DilutedTraceFactory.hpp"
#include "h5-wrapper.hpp"

#include <omp.h>

#include <mutex>

struct DiagramParts {
  DiagramParts(RandomVector const &random_vector,
               Perambulator const &perambulator,
               OperatorFactory const &meson_operator,
               DilutionScheme const &dilution_scheme,
               ssize_t const dilT,
               ssize_t const dilE,
               ssize_t const nev,
               ssize_t const Lt,
               DilutedFactorIndicesCollection const &dil_fac_lookup,
               TraceIndicesCollection const &trace_indices_map)
      : q0(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q0")),
        q1(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q1")),
        q2(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q2")),
        q0q2(q0, q2) {
    trace_factories["trQ1"] = std::unique_ptr<AbstractDilutedTraceFactory>(
        new DilutedTrace1Factory<DilutedFactorType::Q1>(
            q1, trace_indices_map.at("trQ1"), dilution_scheme));

    trace_factories["trQ1Q1"] = std::unique_ptr<AbstractDilutedTraceFactory>(
        new DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1>(
            q1, q1, trace_indices_map.at("trQ1Q1"), dilution_scheme));

    trace_factories["trQ0Q2"] = std::unique_ptr<AbstractDilutedTraceFactory>(
        new DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2>(
            q0, q2, trace_indices_map.at("trQ0Q2"), dilution_scheme));

    trace_factories["trQ1Q1Q1"] = std::unique_ptr<AbstractDilutedTraceFactory>(
        new DilutedTrace3Factory<DilutedFactorType::Q1,
                                 DilutedFactorType::Q1,
                                 DilutedFactorType::Q1>(
            q1, q1, q1, trace_indices_map.at("trQ1Q1Q1"), dilution_scheme));

    trace_factories["trQ1Q0Q2"] = std::unique_ptr<AbstractDilutedTraceFactory>(
        new DilutedTrace3Factory<DilutedFactorType::Q1,
                                 DilutedFactorType::Q0,
                                 DilutedFactorType::Q2>(
            q1, q0, q2, trace_indices_map.at("trQ1Q0Q2"), dilution_scheme));

    trace_factories["trQ1Q1Q1Q1"] = std::unique_ptr<AbstractDilutedTraceFactory>(
        new DilutedTrace4Factory<DilutedFactorType::Q1,
                                 DilutedFactorType::Q1,
                                 DilutedFactorType::Q1,
                                 DilutedFactorType::Q1>(
            q1, q1, q1, q1, q0q2, trace_indices_map.at("trQ1Q1Q1Q1"), dilution_scheme));

    trace_factories["trQ2Q0Q2Q0"] = std::unique_ptr<AbstractDilutedTraceFactory>(
        new DilutedTrace4Factory<DilutedFactorType::Q2,
                                 DilutedFactorType::Q0,
                                 DilutedFactorType::Q2,
                                 DilutedFactorType::Q0>(
            q2, q0, q2, q0, q0q2, trace_indices_map.at("trQ2Q0Q2Q0"), dilution_scheme));

    trace_factories["trQ2Q0Q2Q0Q2Q0"] = std::unique_ptr<AbstractDilutedTraceFactory>(
        new DilutedTrace6Factory<DilutedFactorType::Q2,
                                 DilutedFactorType::Q0,
                                 DilutedFactorType::Q2,
                                 DilutedFactorType::Q0,
                                 DilutedFactorType::Q2,
                                 DilutedFactorType::Q0>(
            q2,
            q0,
            q2,
            q0,
            q2,
            q0,
            q0q2,
            trace_indices_map.at("trQ2Q0Q2Q0Q2Q0"),
            dilution_scheme));
  }

  void clear() {
    q0.clear();
    q1.clear();
    q2.clear();
    q0q2.clear();

    for (auto const &elem : trace_factories) {
      elem.second->clear();
    }
  }

  DilutedFactorFactory<DilutedFactorType::Q0> q0;
  DilutedFactorFactory<DilutedFactorType::Q1> q1;
  DilutedFactorFactory<DilutedFactorType::Q2> q2;

  DilutedProductFactoryQ0Q2 q0q2;

  std::map<std::string, std::unique_ptr<AbstractDilutedTraceFactory>> trace_factories;
};

class Diagram {
 public:
  using AccumulatorVector = std::vector<Accumulator<Complex>>;

  Diagram(std::vector<CorrelatorRequest> const &correlator_requests,
          std::string const &output_path,
          std::string const &output_filename,
          int const Lt,
          std::string const &name)
      : correlator_requests_(correlator_requests),
        output_path_(output_path),
        output_filename_(output_filename),
        Lt_(Lt),
        correlator_(
            Lt, AccumulatorVector(correlator_requests.size(), Accumulator<Complex>{})),
        c_(omp_get_max_threads(),
           AccumulatorVector(correlator_requests.size(), Accumulator<Complex>{})),
        mutexes_(Lt),
        name_(name) {
    assert(output_path_ != "");
    assert(output_filename_ != "");
  }

  std::vector<CorrelatorRequest> const &correlator_requests() const {
    return correlator_requests_;
  }

  void assemble(int const t, BlockIterator const &slice_pair, DiagramParts &q);

  void write();

  std::string const &name() const { return name_; }

  void assemble_impl(AccumulatorVector &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q);

  std::vector<CorrelatorRequest> const &correlator_requests_;

 private:
  std::string const &output_path_;
  std::string const &output_filename_;

  int const Lt_;

  /** OpenMP-shared correlators, indices are (1) time and (2) correlator id. */
  std::vector<AccumulatorVector> correlator_;

  /** OpenMP-shared correlators, indices are (1) thread id and (2) correlator id. */
  std::vector<AccumulatorVector> c_;

  std::vector<std::mutex> mutexes_;

  std::string name_;
};
