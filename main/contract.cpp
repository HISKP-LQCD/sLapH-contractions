/*! @file contract.cpp
 *  Main function of Contraction code for sLapH Perambulators
 */

#include "Correlators.hpp"
#include "git.hpp"
#include "global_data.hpp"
#include "global_data_build_IO_names.hpp"
#include "timings.hpp"

#include <omp.h>

#include <iostream>

/*! Read parameters from infile and perform the specified contractions
 *
 *  In succession instanciate GlobalData, Perambulator,
 *  RandomVector, OperatorsForMesons and Correlators.
 *  - Get paths, physical quantum numbers and desired operators from infile
 *  - Loop over Configuration
 *  - Read perambulators, randomvectors and contract
 *
 *  The flow of this function is depicted in the Flowchart below. The colors
 *  uniquely encode the class a member function belongs to. Beige fields are
 *  not classmembers
 *  @image html contract.pdf "Flow chart of main()"
 */
int main(int ac, char *av[]) {
  TimingScope<1> timing_scope("main");

  std::cout << "This is sLapH-contractions:\n"
            << "  git branch " << git_refspec << "\n"
            << "  git revision " << git_sha1 << "\n"
            << "  git state " << git_changes << "\n"
            << "  compiled by " << git_user << " on " << git_host << "\n"
            << "  running with up to " << omp_get_max_threads() << " OpenMP threads\n";

  // reading global parameters from input file
  GlobalData gd;

  read_parameters(gd, ac, av);

  // initialization of OMP paralization
  Eigen::initParallel();
  Eigen::setNbThreads(gd.nb_eigen_threads);


  OperatorFactory meson_operators(gd.Lt,
                                  gd.Lx,
                                  gd.Ly,
                                  gd.Lz,
                                  gd.number_of_eigen_vec,
                                  (gd.quarks)[0].number_of_dilution_E,
                                  gd.operator_lookuptable,
                                  gd.handling_vdaggerv,
                                  gd.path_vdaggerv,
                                  gd.path_config,
                                  gd.hyp_parameters);

  if (gd.delta_config <= 0) {
    std::cerr
        << "The 'delta_config' option has been set to a non-positive value. This makes "
           "no sense."
        << std::endl;
    std::abort();
  }

  // Loop over all configurations stated in the infile
  for (ssize_t config_i = gd.start_config; config_i <= gd.end_config;
       config_i += gd.delta_config) {
    std::cout << "\nprocessing configuration: " << config_i << "\n\n" << std::endl;
    // changes all paths and names which depend on the configuration
    build_IO_names(gd, config_i);

    // with this odd ordering we can optimize memory usage. Operator creation below
    // will incur a significant temporary memory cost, so we only allocate the bare
    // minimum here (the random vectors)
    RandomVector randomvectors(gd.rnd_vec_construct.nb_entities,
                               gd.rnd_vec_construct.length);

    // read eigenvectors and build operators
    meson_operators.create_operators(gd.filename_eigenvectors, randomvectors, config_i, gd);

    // only now that VdaggerV has been created / loaded do we allocate memory for the
    // perambulators, the allocation overhead is negligible
    Perambulator perambulators(gd.peram_construct.nb_entities,
                               gd.peram_construct.size_rows,
                               gd.peram_construct.size_cols);

    if (gd.handling_vdaggerv != "only_vdaggerv_compute_save"){
      // read perambulators
      perambulators.read_perambulators_from_separate_files(
          gd.Lt, gd.number_of_eigen_vec, gd.quarks, gd.peram_construct.filename_list);
      // read random vectors
      randomvectors.read_random_vectors_from_separate_files(
          gd.rnd_vec_construct.filename_list);
      
      contract(gd.Lt,
              (gd.quarks)[0].number_of_dilution_T,
              (gd.quarks)[0].number_of_dilution_E,
              gd.number_of_eigen_vec,
              meson_operators,
              randomvectors,
              perambulators,
              gd.operator_lookuptable,
              gd.trace_indices_map,
              gd.correlator_requests_map,
              gd.quarkline_lookuptable,
              gd.path_output,
              gd.filename_ending_correlators,
              gd.single_time_slice_combination);
    }
  }
}
