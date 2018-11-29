#include "Correlators.hpp"

//#define DILUTION_ITERATOR_PRINT

#include "Diagram.hpp"
#include "DilutedFactor.hpp"
#include "DilutedFactorY.hpp"
#include "StopWatch.hpp"
#include "dilution-iterator.hpp"
#include "local_timer.hpp"
#include "typedefs.hpp"

#include <omp.h>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

int get_time_delta(BlockIterator const &slice_pair, int const Lt) {
  return abs((slice_pair.sink() - slice_pair.source() - Lt) % Lt);
}

/**
 *  @param quarklines       Instance of Quarklines. Contains prebuilt
 *                          combinations of operators and perambulators
 *  @param meson_operator   Instance of OperatorsForMesons. Contains
 *                          operators (@f$ V^\dagger V $f$) with momenta
 *                          and with/without dilution.
 *  @param perambulators    Instance of Perambulator class. Contains
 *                          Perambulator data
 *  @param operator_lookup
 *  @param corr_lookup
 *  @param quark_lookup
 *
 *  If a diagram is not specified in the infile, corr_lookup contains an empty
 *  vector for this diagram and the build function immediately returns
 */
void contract(const ssize_t Lt,
              const ssize_t dilT,
              const ssize_t dilE,
              const ssize_t nev,
              OperatorFactory const &meson_operator,
              RandomVector const &randomvectors,
              Perambulator const &perambulators,
              OperatorLookup const &operator_lookup,
              DiagramIndicesCollection const &corr_lookup,
              TraceIndicesCollection const &trace_indices_map,
              CorrelatorRequestsMap const &correlator_requests_map,
              DilutedFactorIndicesCollection const &quark_lookup,
              std::string const output_path,
              std::string const output_filename) {
  std::vector<std::unique_ptr<Diagram>> diagrams;

  if (!corr_lookup.at("C2c").empty())
    diagrams.emplace_back(new C2c(corr_lookup.at("C2c"),
                                  correlator_requests_map.at("C2c"),
                                  output_path,
                                  output_filename,
                                  Lt));
  if (!corr_lookup.at("C20").empty())
    diagrams.emplace_back(new C20(corr_lookup.at("C20"),
                                  correlator_requests_map.at("C20"),
                                  output_path,
                                  output_filename,
                                  Lt));
  if (!corr_lookup.at("C20V").empty())
    diagrams.emplace_back(new C20V(corr_lookup.at("C20V"),
                                   correlator_requests_map.at("C20V"),
                                   output_path,
                                   output_filename,
                                   Lt));

  if (!corr_lookup.at("C3c").empty())
    diagrams.emplace_back(new C3c(corr_lookup.at("C3c"),
                                  correlator_requests_map.at("C3c"),
                                  output_path,
                                  output_filename,
                                  Lt));
  if (!corr_lookup.at("C30").empty())
    diagrams.emplace_back(new C30(corr_lookup.at("C30"),
                                  correlator_requests_map.at("C30"),
                                  output_path,
                                  output_filename,
                                  Lt));
  if (!corr_lookup.at("C30V").empty())
    diagrams.emplace_back(new C30V(corr_lookup.at("C30V"),
                                   correlator_requests_map.at("C30V"),
                                   output_path,
                                   output_filename,
                                   Lt));

  if (!corr_lookup.at("C4cB").empty())
    diagrams.emplace_back(new C4cB(corr_lookup.at("C4cB"),
                                   correlator_requests_map.at("C4cB"),
                                   output_path,
                                   output_filename,
                                   Lt));
  if (!corr_lookup.at("C4cC").empty())
    diagrams.emplace_back(new C4cC(corr_lookup.at("C4cC"),
                                   correlator_requests_map.at("C4cC"),
                                   output_path,
                                   output_filename,
                                   Lt));
  if (!corr_lookup.at("C40B").empty())
    diagrams.emplace_back(new C40B(corr_lookup.at("C40B"),
                                   correlator_requests_map.at("C40B"),
                                   output_path,
                                   output_filename,
                                   Lt));
  if (!corr_lookup.at("C40C").empty())
    diagrams.emplace_back(new C40C(corr_lookup.at("C40C"),
                                   correlator_requests_map.at("C40C"),
                                   output_path,
                                   output_filename,
                                   Lt));

  if (!corr_lookup.at("C4cD").empty())
    diagrams.emplace_back(new C4cD(corr_lookup.at("C4cD"),
                                   correlator_requests_map.at("C4cD"),
                                   output_path,
                                   output_filename,
                                   Lt));
  if (!corr_lookup.at("C4cV").empty())
    diagrams.emplace_back(new C4cV(corr_lookup.at("C4cV"),
                                   correlator_requests_map.at("C4cV"),
                                   output_path,
                                   output_filename,
                                   Lt));
  if (!corr_lookup.at("C40D").empty())
    diagrams.emplace_back(new C40D(corr_lookup.at("C40D"),
                                   correlator_requests_map.at("C40D"),
                                   output_path,
                                   output_filename,
                                   Lt));
  if (!corr_lookup.at("C40V").empty())
    diagrams.emplace_back(new C40V(corr_lookup.at("C40V"),
                                   correlator_requests_map.at("C40V"),
                                   output_path,
                                   output_filename,
                                   Lt));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  StopWatch swatch("All contractions");

#pragma omp parallel
  {
    swatch.start();

    LT_CORRELATOR_DECLARE;

    DiagramParts q(randomvectors,
                   perambulators,
                   meson_operator,
                   dilution_scheme,
                   dilT,
                   dilE,
                   nev,
                   Lt,
                   quark_lookup,
                   corr_lookup,
                   trace_indices_map);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
#pragma omp critical(cout)
      {
        std::cout << "Thread " << std::setw(3) << omp_get_thread_num() << " of "
                  << std::setw(3) << omp_get_num_threads() << " starts with block pair "
                  << std::setw(5) << b << " of " << std::setw(5) << dilution_scheme.size()
                  << "." << std::endl;
      }

      auto const block_pair = dilution_scheme[b];

      // Build the diagrams.
      for (auto &diagram : diagrams) {
        if (diagram->corr_lookup().empty()) {
          continue;
        }
        LT_CORRELATOR_START;

        for (auto const slice_pair : block_pair) {
          int const t = get_time_delta(slice_pair, Lt);

          diagram->assemble(t, slice_pair, q);
        }  // End of slice pair loop.
        LT_CORRELATOR_STOP;
        LT_CORRELATOR_PRINT(std::string("[contract] ") + diagram->name());
      }  // End of diagram loop.

      q.clear();
    }  // End of block pair loop.

    swatch.stop();
  }  // End of parallel section.

  swatch.print();

  for (auto &diagram : diagrams) {
    diagram->write();
  }
}
