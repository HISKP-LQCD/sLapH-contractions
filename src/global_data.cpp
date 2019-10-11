#include "global_data.hpp"

#include "git.hpp"
#include "global_data_build_IO_names.hpp"
#include "global_data_input_handling.hpp"
#include "init_lookup_tables.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

namespace po = boost::program_options;

GlobalData::GlobalData() {
  momentum_cutoff[0] = 4;
  momentum_cutoff[1] = 5;
  momentum_cutoff[2] = 6;
  momentum_cutoff[3] = 7;
  momentum_cutoff[4] = 4;

  quarkline_lookuptable["Q0"];
  quarkline_lookuptable["Q1"];
  quarkline_lookuptable["Q2"];

  trace_indices_map["trQ0Q2"];
  trace_indices_map["trQ1"];
  trace_indices_map["trQ1Q0Q2"];
  trace_indices_map["trQ1Q1"];
  trace_indices_map["trQ1Q1Q1"];
  trace_indices_map["trQ1Q1Q1Q1"];
  trace_indices_map["trQ2Q0Q2Q0"];
  trace_indices_map["trQ2Q0Q2Q0Q2Q0"];
}

/**
 * Convenience function for when a \ref store_to value is being provided to
 * typed_value.
 *
 * @param store_to The variable that will hold the parsed value upon notify.
 *
 * @return Pointer to a type_value.
 */
template <typename T>
boost::program_options::typed_value<T> *make_value(T *store_to) {
  return boost::program_options::value<T>(store_to);
}

void read_parameters(GlobalData &gd, int ac, char *av[]) {
  std::string input_file;
  std::string output_file;

  // Variables that will store parsed values for quarks.
  std::vector<std::string> quark_configs;
  // Variables that will store parsed values for operators.
  std::vector<std::string> operator_list_configs;
  // Variables that will store parsed values for correlators.
  std::vector<std::string> correlator_list_configs;

  // Declare a group of options that will be allowed only on command line
  po::options_description generic("Command line options");
  generic.add_options()("help,h", "produce help message");
  generic.add_options()("version,v", "print version string");
  generic.add_options()("verbose", "does additional tests and prints more details");
  generic.add_options()("input,i",
                        po::value<std::string>(&input_file)->default_value("LapHs.in"),
                        "name of input file.");
  generic.add_options()("output,o",
                        po::value<std::string>(&output_file)->default_value("LapHs.out"),
                        "name of output file.");

  generic.add_options()(
      "single_time_slice_combination",
      po::value<int>(&gd.single_time_slice_combination)->default_value(-1),
      "Only compute this time slice combination (for benchmarking only)");

  // Declare a group of options that will be allowed both on command line and
  // in input file
  po::options_description config("Input file options");

  //////////////////////////////////////////////////////////////////////////////
  // Options for infile ////////////////////////////////////////////////////////

  // parallelisation options
  config.add_options()("nb_eigen_threads",
                       po::value<ssize_t>(&gd.nb_eigen_threads)->default_value(1),
                       "nb_eigen_threads: number of threads Eigen uses internally");

  // lattice options
  config.add_options()(
      "output_path",
      po::value<std::string>(&gd.path_output)->default_value("../../contractions"),
      "path for output");
  config.add_options()("path_config",
                       po::value<std::string>(&gd.path_config),
                       "path for gauge configurations");
  config.add_options()(
      "lattice", po::value<std::string>(&gd.name_lattice), "Codename of the lattice");
  config.add_options()(
      "Lt", po::value<int>(&gd.Lt)->default_value(0), "Lt: temporal lattice extent");
  config.add_options()("Lx",
                       po::value<int>(&gd.Lx)->default_value(0),
                       "Lx: lattice extent in x direction");
  config.add_options()("Ly",
                       po::value<int>(&gd.Ly)->default_value(0),
                       "Ly: lattice extent in y direction");
  config.add_options()("Lz",
                       po::value<int>(&gd.Lz)->default_value(0),
                       "Lz: lattice extent in z direction");

  config.add_options()("alpha1",
                       po::value<double>(&gd.hyp_parameters.alpha1)->default_value(0),
                       "alpha1: Inner level strength of Hyp-smearing");
  config.add_options()("alpha2",
                       po::value<double>(&gd.hyp_parameters.alpha2)->default_value(0),
                       "alpha2: Outer level strength of Hyp-smearing");
  config.add_options()(
      "iterations",
      po::value<ssize_t>(&gd.hyp_parameters.iterations)->default_value(0),
      "iterations: Number of Hyp-smearing applications");

  // eigenvector options
  config.add_options()("number_of_eigen_vec",
                       po::value<int>(&gd.number_of_eigen_vec)->default_value(0),
                       "Number of eigen vectors");
  config.add_options()("path_eigenvectors",
                       po::value<std::string>(&gd.path_eigenvectors)->default_value("."),
                       "directory of eigenvectors");
  config.add_options()(
      "name_eigenvectors",
      po::value<std::string>(&gd.name_eigenvectors)->default_value("eigenvector"),
      "name of eigenvectors\nThe full name is internally created to:\n"
      "\"name_of_eigenvectors.configuration.time_slice\"");
  config.add_options()(
      "handling_vdaggerv",
      po::value<std::string>(&gd.handling_vdaggerv)->default_value("build"),
      "The options are:\n"
      "build: VdaggerV is build for all operators but not written to disk\n"
      "write: VdaggerV is build for all operators and written to disk\n"
      "read: VdaggerV was previously constructed and is read from disk");
  config.add_options()("path_vdaggerv",
                       po::value<std::string>(&gd.path_vdaggerv)->default_value(""),
                       "Path of vdaggerv");

  // quark options
  config.add_options()("quarks.quark",
                       make_value(&quark_configs),
                       "quark input, must be of type:\n"
                       "quark = \n type:number of rnd. vec.:\n"
                       " dil type time:number of dil time:\n"
                       " dil type ev:number of dil ev:\n"
                       " dil type Dirac:number of dil Dirac:\n"
                       " path of the perambulators for these quarks");

  // operator list options
  config.add_options()("operator_lists.operator_list",
                       make_value(&operator_list_configs),
                       "operator input is rather complicated - see documentation!!");

  // correlator list options
  config.add_options()("correlator_lists.correlator_list",
                       make_value(&correlator_list_configs),
                       "correlator input is rather complicated - see documentation!!");

  // configuration options
  config.add_options()("start_config",
                       po::value<int>(&gd.start_config)->default_value(-1),
                       "First configuration");
  config.add_options()("end_config",
                       po::value<int>(&gd.end_config)->default_value(0),
                       "Last configuration");
  config.add_options()("delta_config",
                       po::value<int>(&gd.delta_config)->default_value(1),
                       "Stepsize between two configurations");

  config.add_options()("momentum_cutoff_0",
                       po::value<int>(&gd.momentum_cutoff[0])->default_value(gd.momentum_cutoff[0]),
                       "Cutoff for \\sum_i |pi|² when |P|² = 0");
  config.add_options()("momentum_cutoff_1",
                       po::value<int>(&gd.momentum_cutoff[1])->default_value(gd.momentum_cutoff[1]),
                       "Cutoff for \\sum_i |pi|² when |P|² = 1");
  config.add_options()("momentum_cutoff_2",
                       po::value<int>(&gd.momentum_cutoff[2])->default_value(gd.momentum_cutoff[2]),
                       "Cutoff for \\sum_i |pi|² when |P|² = 2");
  config.add_options()("momentum_cutoff_3",
                       po::value<int>(&gd.momentum_cutoff[3])->default_value(gd.momentum_cutoff[3]),
                       "Cutoff for \\sum_i |pi|² when |P|² = 3");
  config.add_options()("momentum_cutoff_4",
                       po::value<int>(&gd.momentum_cutoff[4])->default_value(gd.momentum_cutoff[4]),
                       "Cutoff for \\sum_i |pi|² when |P|² = 4");

  //////////////////////////////////////////////////////////////////////////////

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config);

  po::options_description input_file_options;
  input_file_options.add(config);

  po::options_description visible("Allowed options");
  visible.add(generic).add(config);
  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(),
            vm);
  po::notify(vm);

  // command line options
  if (vm.count("help")) {
    std::cout << visible << "\n";
    exit(0);
  }
  if (vm.count("verbose")) {
    gd.verbose = 1;
  } else
    gd.verbose = 0;
  if (vm.count("version")) {
    std::cout << "Contraction code for LapHs perambulators.\n"
              << "Git SHA-1: " << git_sha1 << "\n"
              << "Git refspec: " << git_refspec << "\n"
              << "Git status: " << git_changes << "\n";
    exit(0);
  }

  {
    std::cout << "\nConfiguration file read in is:\n\n";

    std::ifstream ifs(input_file.c_str());
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    std::cout << buffer.str() << "\n\n---- End of configuration file ----\n\n";
  }

  {
    std::ifstream ifs(input_file.c_str());
    if (!ifs) {
      std::cout << "CANNOT open input file: " << input_file << "\n";
      exit(1);
    } else {
      po::store(parse_config_file(ifs, input_file_options), vm);
      po::notify(vm);
    }
  }

  /****************************************************************************/

  // checks, terminal output and munging of strings for quarks, operators and
  // correlators
  input_handling(gd, quark_configs, operator_list_configs, correlator_list_configs);

  // setting the lookup tables for all needed quantum numbers to calculate
  // the wanted correlators
  init_lookup_tables(gd);

  std::cout << "The following elements are defined in `gd.trace_indices_map`:\n";
  for (auto const &elem : gd.trace_indices_map)
    std::cout << "  - " << elem.first << "\n";
  std::cout << std::endl;

  // setting the sizes and numbers of random vectors and perambulators
  /**
   * @todo: setting the sizes and numbers of rnd_vecs and perams should be put
   * in a separate function
   */
  gd.rnd_vec_construct.nb_entities = 0;
  for (const auto &q : gd.quarks)
    gd.rnd_vec_construct.nb_entities += q.number_of_rnd_vec;
  gd.rnd_vec_construct.length = gd.Lt * 4 * gd.number_of_eigen_vec;
  gd.peram_construct.nb_entities = gd.rnd_vec_construct.nb_entities;
  for (const auto &q : gd.quarks) {
    for (ssize_t r = 0; r < q.number_of_rnd_vec; r++) {
      gd.peram_construct.size_rows.push_back(gd.rnd_vec_construct.length);
      gd.peram_construct.size_cols.push_back((gd.Lt / q.number_of_dilution_T) *
                                             q.number_of_dilution_E *
                                             q.number_of_dilution_D);
    }
  }

  std::cout << gd << std::endl;

  // printing information about memory consumption of all relevant parts that are cached
  std::cout << "Memory consumption:" << std::endl;

  std::cout << "\tOperatorFactory:\t" << std::fixed << std::setprecision(2)
            << gd.operator_lookuptable.size() * gd.Lt * gd.number_of_eigen_vec *
                   gd.number_of_eigen_vec * sizeof(Complex) / std::pow(2, 30)
            << " Gb" << std::endl;

  std::cout << "\tRandomVector:\t" << std::fixed << std::setprecision(2)
            << gd.rnd_vec_construct.nb_entities * gd.rnd_vec_construct.length *
                   sizeof(Complex) / std::pow(2, 30)
            << " Gb" << std::endl;

  int peram_matrix_size_sum = 0;
  for (auto i = 0; i < ssize(gd.peram_construct.size_rows); ++i) {
    peram_matrix_size_sum +=
        gd.peram_construct.size_rows[i] * gd.peram_construct.size_cols[i];
  }
  std::cout << "\tPerambulator:\t" << std::fixed << std::setprecision(2)
            << peram_matrix_size_sum * sizeof(Complex) / std::pow(2, 30) << " Gb"
            << std::endl;

  std::cout << "\tDiagramParts:" << std::endl;

  int total_number_of_random_combinations_in_Q0 = 0;
  for (auto const &q : gd.quarkline_lookuptable.at("Q0")) {
    total_number_of_random_combinations_in_Q0 += q.rnd_vec_ids.size();
  }
  int Q0_matrix_size =
      gd.quarks[0].number_of_dilution_D * gd.quarks[0].number_of_dilution_E;
  std::cout << "\t\tQ0:\t" << std::fixed << std::setprecision(2)
            << gd.Lt / gd.quarks[0].number_of_dilution_T *
                   (gd.Lt / gd.quarks[0].number_of_dilution_T - 1) / 2 *
                   total_number_of_random_combinations_in_Q0 * Q0_matrix_size *
                   sizeof(Complex) / std::pow(2, 30)
            << " Gb" << std::endl;

  int total_number_of_random_combinations_in_Q1 = 0;
  for (auto const &q : gd.quarkline_lookuptable.at("Q1")) {
    total_number_of_random_combinations_in_Q1 += q.rnd_vec_ids.size();
  }
  int Q1_matrix_size =
      gd.quarks[0].number_of_dilution_D * gd.quarks[0].number_of_dilution_E;
  std::cout << "\t\tQ1:\t" << std::fixed << std::setprecision(2)
            << gd.Lt / gd.quarks[0].number_of_dilution_T *
                   (gd.Lt / gd.quarks[0].number_of_dilution_T - 1) / 2 *
                   total_number_of_random_combinations_in_Q1 * Q1_matrix_size *
                   sizeof(Complex) / std::pow(2, 30)
            << " Gb" << std::endl;

  int total_number_of_random_combinations_in_Q2 = 0;
  for (auto const &q : gd.quarkline_lookuptable.at("Q2")) {
    total_number_of_random_combinations_in_Q2 += q.rnd_vec_ids.size();
  }
  int Q2_matrix_size =
      gd.quarks[0].number_of_dilution_D * gd.quarks[0].number_of_dilution_E;
  std::cout << "\t\tQ2:\t" << std::fixed << std::setprecision(2)
            << gd.Lt / gd.quarks[0].number_of_dilution_T *
                   (gd.Lt / gd.quarks[0].number_of_dilution_T - 1) / 2 *
                   total_number_of_random_combinations_in_Q2 * Q2_matrix_size *
                   sizeof(Complex) / std::pow(2, 30)
            << " Gb" << std::endl;

  int total_number_of_random_combinations_in_trQ1Q1 = 0;
  for (auto const &q : gd.trace_indices_map.at("trQ1Q1")) {
    total_number_of_random_combinations_in_trQ1Q1 +=
        gd.quarkline_lookuptable.at("Q1").at(q[0]).rnd_vec_ids.size() *
        gd.quarkline_lookuptable.at("Q1").at(q[1]).rnd_vec_ids.size();
  }
  std::cout << "\ttrQ1Q1:\t" << std::fixed << std::setprecision(2)
            << gd.Lt * gd.Lt * total_number_of_random_combinations_in_trQ1Q1 *
                   sizeof(Complex) / std::pow(2, 30)
            << " Gb" << std::endl;

  int total_number_of_random_combinations_in_trQ0Q2 = 0;
  for (auto const &q : gd.trace_indices_map.at("trQ0Q2")) {
    auto const size1 = gd.quarkline_lookuptable.at("Q0").at(q[0]).rnd_vec_ids.size();
    auto const size2 = gd.quarkline_lookuptable.at("Q2").at(q[1]).rnd_vec_ids.size();
    total_number_of_random_combinations_in_trQ0Q2 += size1 * size2;
  }
  std::cout << "\ttrQ0Q2:\t" << std::fixed << std::setprecision(2)
            << gd.Lt * gd.Lt * total_number_of_random_combinations_in_trQ0Q2 *
                   sizeof(Complex) / std::pow(2, 30)
            << " Gb" << std::endl;

  int total_number_of_random_combinations_in_trQ1 = 0;
  for (auto const &q : gd.trace_indices_map.at("trQ1")) {
    total_number_of_random_combinations_in_trQ1 +=
        gd.quarkline_lookuptable.at("Q1").at(q[0]).rnd_vec_ids.size();
  }
  std::cout << "\ttrQ1:\t" << std::fixed << std::setprecision(2)
            << gd.Lt * total_number_of_random_combinations_in_trQ1 * sizeof(Complex) /
                   std::pow(2, 30)
            << " Gb" << std::endl;

  std::cout << "\tDiagrams:" << std::endl;
}

#define GLOBAL_DATA_PRINT(x) (std::cout << "    " << #x << ": " << x << "\n")

std::ostream &operator<<(std::ostream &os, std::map<int, int> const &map) {
  for (auto const &elem : map) {
    os << elem.first << "->" << elem.second << " ";
  }
  return os;
}

std::ostream &operator<<(std::ostream &os, GlobalData const &gd) {
  std::cout << "\nGlobal Data:\n";

  std::cout << "  Lattice size:\n";
  GLOBAL_DATA_PRINT(gd.Lx);
  GLOBAL_DATA_PRINT(gd.Ly);
  GLOBAL_DATA_PRINT(gd.Lz);
  GLOBAL_DATA_PRINT(gd.Lt);

  std::cout << "  Other stuff:\n";
  GLOBAL_DATA_PRINT(gd.dim_row);
  GLOBAL_DATA_PRINT(gd.V_TS);
  GLOBAL_DATA_PRINT(gd.V_for_lime);
  GLOBAL_DATA_PRINT(gd.number_of_eigen_vec);
  GLOBAL_DATA_PRINT(gd.number_of_inversions);
  GLOBAL_DATA_PRINT(gd.start_config);
  GLOBAL_DATA_PRINT(gd.end_config);
  GLOBAL_DATA_PRINT(gd.delta_config);
  GLOBAL_DATA_PRINT(gd.verbose);
  GLOBAL_DATA_PRINT(gd.nb_eigen_threads);
  GLOBAL_DATA_PRINT(gd.path_eigenvectors);
  GLOBAL_DATA_PRINT(gd.name_eigenvectors);
  GLOBAL_DATA_PRINT(gd.filename_eigenvectors);
  GLOBAL_DATA_PRINT(gd.path_perambulators);
  GLOBAL_DATA_PRINT(gd.name_perambulators);
  GLOBAL_DATA_PRINT(gd.name_lattice);
  GLOBAL_DATA_PRINT(gd.filename_ending_correlators);
  GLOBAL_DATA_PRINT(gd.path_output);
  GLOBAL_DATA_PRINT(gd.path_config);
  GLOBAL_DATA_PRINT(gd.handling_vdaggerv);
  GLOBAL_DATA_PRINT(gd.path_vdaggerv);

  std::cout << "  Momentum Cutoff: ";
  std::cout << gd.momentum_cutoff << "\n";

  //! @TODO Print more stuff here.

  return os;
}
