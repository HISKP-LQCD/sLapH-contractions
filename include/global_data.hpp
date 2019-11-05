#pragma once

#include "global_data_typedefs.hpp"
#include "typedefs.hpp"

#include <sys/stat.h>
#include <array>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

typedef std::map<std::string, std::vector<DiagramIndex>> DiagramIndicesCollection;
typedef std::map<std::string, std::vector<Indices>> TraceIndicesCollection;

enum class Location { source, sink };

bool operator<(std::vector<Location> const &left, std::vector<Location> const &right);

struct TraceRequest {
  std::string tr_name;
  ssize_t tr_id;
  std::vector<Location> locations;

  bool operator==(TraceRequest const &other) const {
    return tr_name == other.tr_name && tr_id == other.tr_id &&
           locations == other.locations;
  }

  bool operator<(TraceRequest const &other) const {
    return tr_name < other.tr_name && tr_id < other.tr_id && locations < other.locations;
  }
};

struct CorrelatorRequest {
  std::vector<TraceRequest> trace_requests;
  std::string hdf5_dataset_name;

  bool operator==(CorrelatorRequest const &other) const {
    return trace_requests == other.trace_requests;
  }
};

typedef std::map<std::string, std::vector<CorrelatorRequest>> CorrelatorRequestsMap;

/**
 * Class containing all metadata for contractions and functions to set them
 * from infile
 *
 *  Metadata roughly characterized by either
 *
 *  - physical parameters
 *  - flags
 *  - paths
 */
struct GlobalData {
  GlobalData();

  int Lx, Ly, Lz, Lt;
  int dim_row, V_TS, V_for_lime;
  int number_of_eigen_vec;
  int number_of_inversions;
  int start_config, end_config, delta_config;
  int verbose;
  int max_momentum;
  ssize_t nb_eigen_threads;

  ssize_t nb_evec_read_threads;
  ssize_t nb_vdaggerv_eigen_threads;

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
  std::string path_correlator_list;

  RandomVectorConstruction rnd_vec_construct;
  PerambulatorConstruction peram_construct;

  std::vector<quark> quarks;
  Operator_list operator_list;
  Correlator_list correlator_list;
  DilutedFactorIndicesCollection quarkline_lookuptable;
  OperatorLookup operator_lookuptable;
  TraceIndicesCollection trace_indices_map;
  CorrelatorRequestsMap correlator_requests_map;

  HypPars hyp_parameters;

  int single_time_slice_combination;
};

/**
 * Reading the input parameters from the infile in the main routine and
 * initializing GlobalData.
 */
void read_parameters(GlobalData &gd, int ac, char *av[]);

std::ostream &operator<<(std::ostream &os, GlobalData const &gd);
