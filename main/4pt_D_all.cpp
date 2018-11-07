// *****************************************************************************
// Name         : contract.cpp
// Author       : Bastian Knippschild
// Contributors : Christopher Helmes, Christian Jost, Liuming Liu, Markus Werner
// Version      : 0.1
// Copyright    : Copies are prohibited so far
// Description  : Contraction code for Laphs Perambulators
// *****************************************************************************

#include "global_data_typedefs.hpp"
#include "Perambulator.hpp"
#include "RandomVector.hpp"
#include "typedefs.hpp"

#include <omp.h>
#include <boost/format.hpp>

#include <iostream>

int main (int ac, char* av[]) {

  // ---------------------------------------------------------------------------
  // Some variables definitions which should be read from infile!
  const int Lt = 48;
  const int Lx = 24;
  const int Ly = 24;
  const int Lz = 24;

  const ssize_t nb_ev = 120;

  const ssize_t nb_peram = 6;
  const ssize_t nb_rnd_vec = 6;
  const ssize_t nb_dil_T = 2;
  const ssize_t nb_dil_E = 6;
  const std::string peram_path = 
                    "/data/LapHs/contraction_Markus/test_data/peram/up";
  const quark q("u", nb_peram, "B", nb_dil_T, "I", nb_dil_E, "F", 4, 0, 
                peram_path);

  const ssize_t start_config = 714;
  const ssize_t end_config = 714;
  const ssize_t delta_config = 1;

  // ---------------------------------------------------------------------------
  // initialization of OMP paralization
  const ssize_t nb_omp_threads = 4;
  const ssize_t nb_eigen_threads = 1;
  Eigen::initParallel();
  omp_set_dynamic(0);
  omp_set_num_threads(nb_omp_threads);
  Eigen::setNbThreads(nb_eigen_threads);

  // ---------------------------------------------------------------------------
  // Creating instances of perambulators, random vectors, operators, and 
  // correlators. The eigenvectors are read from disc in the operator class.
  Perambulator perambulators(nb_peram, 
                     std::vector<ssize_t>(nb_rnd_vec, Lt*4*nb_ev), 
                     std::vector<ssize_t>(nb_rnd_vec, (Lt/nb_dil_T)*nb_dil_E*4));
  RandomVector randomvectors(nb_peram, Lt*4*nb_ev);

  // ---------------------------------------------------------------------------
  // Loop over all configurations stated in the infile -------------------------
  for(ssize_t config_i  = start_config; config_i <= end_config; 
                                       config_i += delta_config){
    std::cout << "\nprocessing configuration: " << config_i << "\n\n";

    // creating names for perambulators and random vectors ---------------------
    std::vector<std::string> peram_names;
    std::vector<std::string> rnd_names;

    for(int rnd_vec_i = 0; rnd_vec_i < nb_rnd_vec; ++rnd_vec_i){
      char temp1[200];
      // building paths and filenames for rnd vecs
      sprintf(temp1, "cnfg%zu/rnd_vec_%01d/", config_i, rnd_vec_i);
      char temp2[200];
      sprintf(temp2, "randomvector.rndvecnb%02d.u.nbev%04zu.%04zu", 
                     rnd_vec_i, nb_ev, config_i);
      rnd_names.push_back(q.path + "/" + temp1 + temp2);
      char temp3[200];
      sprintf(temp3,
          "perambulator.rndvecnb%02d.u.TsoB%04d.VsoI%04zu.DsoF%1d.TsiF%04d."
          "SsiF%d.DsiF4.CsiF3.smeared1.%05zu", 
          rnd_vec_i, int(Lt / nb_dil_T), nb_dil_E,
          4, Lt, Lx*Ly*Lz, config_i);
      peram_names.push_back(q.path + "/" + temp1 + temp3);
    }
    // read perambulators ------------------------------------------------------
    perambulators.read_perambulators_from_separate_files(Lt, nb_ev, 
                                          std::vector<quark>(1,q), peram_names);
    // read random vectors -----------------------------------------------------
    randomvectors.read_random_vectors_from_separate_files(rnd_names);

    // creating one quarkline Q2V and closing it with all other operators!


  }
  // That's all Folks!
  return 0;
}

