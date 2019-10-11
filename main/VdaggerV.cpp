#include "Correlators.hpp"
#include "git.hpp"
#include "global_data.hpp"
#include "global_data_build_IO_names.hpp"
#include "timings.hpp"

#include <omp.h>

#include <iostream>
#include "OperatorsForMesons.hpp"

/*! @file vdaggerv.cpp
 *  *  Main function of creating vdaggerv objects for sLapH contractions
 */
int main(int ac, char *av[]) {
  // Some variables definitions which should be read from infile!
  //
  std::cout << "This is sLapH-contractions:\n"
            << "  git branch " << git_refspec << "\n"
            << "  git revision " << git_sha1 << "\n"
            << "  git state " << git_changes << "\n"
            << "  compiled by " << git_user << " on " << git_host << "\n"
            << "  running with up to " << omp_get_max_threads() << " OpenMP threads\n";

  // reading global parameters from input file
  GlobalData gd;
  read_parameters(gd, ac, av);
  // std::cout<<"Parameters are read"<<std::endl;

  // initialization of OMP paralization
  Eigen::initParallel();
  Eigen::setNbThreads(gd.nb_eigen_threads);

  // Creating lookuptable for Operator construction by hand. TODO: Include
  // this in input file handling!

  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup;

  std::vector<std::pair<char, char>> displacement;
  ssize_t id = 0;
  for (int p1 = -gd.max_momentum; p1 <= gd.max_momentum; p1++)
    for (int p2 = -gd.max_momentum; p2 <= gd.max_momentum; p2++)
      for (int p3 = -gd.max_momentum; p3 <= gd.max_momentum; p3++) {
        int maxQsq = gd.max_momentum * gd.max_momentum;
        if (p1 * p1 + p2 * p2 + p3 * p3 <= maxQsq && !(p1 == 0 && p2 == 0 && p3 == 0)) {
          std::array<int, 3> momentum = {p1, p2, p3};
          auto it = std::find_if(
              gd.operator_lookuptable.vdaggerv_lookup.begin(),
              gd.operator_lookuptable.vdaggerv_lookup.end(),
              [&momentum](VdaggerVQuantumNumbers vdv_qn) {
                const std::array<int, 3> pm = {-momentum[0], -momentum[1], -momentum[2]};
                if (vdv_qn.momentum == pm)
                  return true;
                else
                  return false;
              });
          if (!(it != gd.operator_lookuptable.vdaggerv_lookup.end())) {
            gd.operator_lookuptable.vdaggerv_lookup.emplace_back(VdaggerVQuantumNumbers(
                id, {momentum[0], momentum[1], momentum[2]}, displacement));
            id++;
          }
        }
      }

  // Loop over all configurations stated in the infile -------------------------
  for (ssize_t config_i = gd.start_config; config_i <= gd.end_config;
       config_i += gd.delta_config) {
    std::cout << "\nprocessing configuration: " << config_i << "\n\n";

    build_IO_names(gd, config_i);

    OperatorFactory meson_operators(
        gd.Lt,
        gd.Lx,
        gd.Ly,
        gd.Lz,
        gd.number_of_eigen_vec,
        0,  // We set the eigenvector dilution argument to zero
        gd.operator_lookuptable,
        gd.handling_vdaggerv,
        gd.path_vdaggerv,
        gd.path_config,
        gd.hyp_parameters);

    //  gd.operator_lookuptable=operator_lookuptable;
    meson_operators.build_vdaggerv(gd.filename_eigenvectors, config_i, gd);
  }
}
