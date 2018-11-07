/*! @file global_data.h
 */

#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>
#include <iosfwd>

#include <sys/stat.h>

#include "global_data_typedefs.h"
#include "typedefs.h"

/*! Class containing all metadata for contractions and functions to set them
 *  from infile
 *
 *  Metadata roughly characterized by either
 *  - physical parameters
 *  - flags
 *  - paths
 */
struct GlobalData {
  GlobalData() {
    momentum_cutoff[0] = 4;
    momentum_cutoff[1] = 5;
    momentum_cutoff[2] = 6;
    momentum_cutoff[3] = 7;
    momentum_cutoff[4] = 4;
  }

  int Lx, Ly, Lz, Lt;
  int dim_row, V_TS, V_for_lime;
  int number_of_eigen_vec;
  int number_of_inversions;
  int start_config, end_config, delta_config;
  int verbose;
  ssize_t nb_eigen_threads;
  std::string path_eigenvectors;
  std::string name_eigenvectors;
  std::string filename_eigenvectors;
  std::string path_perambulators;
  std::string name_perambulators;
  std::string name_lattice;
  std::string filename_ending_correlators;
  std::string path_output;
  std::string path_config;
  std::string handling_vdaggerv;
  std::string path_vdaggerv;

  RandomVectorConstruction rnd_vec_construct;
  PerambulatorConstruction peram_construct;

  std::vector<quark> quarks;
  Operator_list operator_list;
  Correlator_list correlator_list;
  DilutedFactorIndicesCollection quarkline_lookuptable;
  OperatorLookup operator_lookuptable;
  DiagramIndicesCollection correlator_lookuptable;

  std::map<int, int> momentum_cutoff;

  HypPars hyp_parameters;
};

/*! Reading the input parameters from the infile in the main routine and
 *  initializing GlobalData
 */
void read_parameters(GlobalData &gd, int ac, char *av[]);

/*! Fills the random vector and perambulator structs with the paths and
 *  file names to read the data
 */
void build_IO_names(GlobalData &gd, const ssize_t config);
  

void input_handling(GlobalData &gd,
                    const std::vector<std::string> &quark_configs,
                    const std::vector<std::string> &operator_list_configs,
                    const std::vector<std::string> &correlator_list_configs);

void init_lookup_tables(GlobalData &gd);


std::ostream &operator<<(std::ostream &os, GlobalData const &gd);
