#include "Correlators.hpp"

#include "Diagram.hpp"
#include "DilutedFactor.hpp"
#include "DilutedFactorFactory.hpp"
#include "StopWatch.hpp"
#include "dilution-iterator.hpp"
#include "timings.hpp"
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
              TraceIndicesCollection const &trace_indices_map,
              CorrelatorRequestsMap const &correlator_requests_map,
              DilutedFactorIndicesCollection const &quark_lookup,
              std::string const output_path,
              std::string const output_filename,
              int single_time_slice_combination) {
  TimingScope<1> timing_scope("contract");

  std::vector<Diagram> diagrams;

  for (auto const &elem : correlator_requests_map) {
    auto const name = elem.first;
    if (!elem.second.empty()) {
      diagrams.emplace_back(elem.second, output_path, output_filename, Lt, name);
    }
  }

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  StopWatch swatch("All contractions");

  swatch.start();

  DiagramParts q(randomvectors,
                 perambulators,
                 meson_operator,
                 dilution_scheme,
                 dilT,
                 dilE,
                 nev,
                 Lt,
                 quark_lookup,
                 trace_indices_map);

  for (int b = 0; b < dilution_scheme.size(); ++b) {
    if (single_time_slice_combination >= 0 && single_time_slice_combination != b) {
      continue;
    }

    std::cout << "Starts with block pair " << std::setw(5) << b << " of " << std::setw(5)
              << dilution_scheme.size() << "." << std::endl;

    auto const block_pair = dilution_scheme[b];

    // Build the diagrams.
    for (auto &diagram : diagrams) {
      if (diagram.correlator_requests().empty()) {
        continue;
      }
      TimingScope<1> timing_scope("contract diagram", diagram.name());

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        diagram.assemble(t, slice_pair, q);
      }  // End of slice pair loop.
    }    // End of diagram loop.

    q.clear();
  }  // End of block pair loop.

  swatch.stop();
  swatch.print();

  for (auto &diagram : diagrams) {
    diagram.write();
  }
}
