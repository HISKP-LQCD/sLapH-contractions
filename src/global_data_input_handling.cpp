#include "global_data_input_handling.hpp"

#include "global_data_utils.hpp"

namespace gdu = ::global_data_utils;

namespace {

using gdu::make_correlator;
using gdu::make_operator_list;
using gdu::make_quark;
using gdu::quark_check;

/**
 * Simplifies and cleans GlobalData::read_parameters().
 *
 * @todo Does this actually make it easier?
 */
void lattice_input_data_handling(const std::string path_output,
                                 const std::string name_lattice,
                                 const std::string path_config,
                                 double alpha1,
                                 double alpha2,
                                 size_t iterations,
                                 int Lt,
                                 int Lx,
                                 int Ly,
                                 int Lz) {
  try {
    if (Lt < 1) {
      std::cout << "\ninput file error:\n"
                << "\toption \"Lt\""
                << " is mendatory and its value must be an integer greater than 0!"
                << "\n\n";
      exit(1);
    } else
      std::cout << "\n\ttemporal lattice extent .................. " << Lt << "\n";
    //
    if (Lx < 1) {
      std::cout << "\ninput file error:\n"
                << "\toption \"Lx\""
                << " is mandatory and its value must be an integer greater than 0!"
                << "\n\n";
      exit(1);
    } else
      std::cout << "\tspatial lattice extent in x direction .... " << Lx << "\n";
    //
    if (Ly < 1) {
      std::cout << "\ninput file error:\n"
                << "\toption \"Ly\""
                << " is mandatory and its value must be an integer greater than 0!"
                << "\n\n";
      exit(1);
    } else
      std::cout << "\tspatial lattice extent in y direction .... " << Ly << "\n";
    //
    if (Lz < 1) {
      std::cout << "\ninput file error:\n"
                << "\toption \"Lz\""
                << " is mandatory and its value must be an integer greater than 0!"
                << "\n\n\n";
      exit(1);
    } else
      std::cout << "\tspatial lattice extent in z direction .... " << Lz << "\n\n";
    std::cout << "\tEnsemble ...................................... " << name_lattice
              << std::endl;
    std::cout << "\tResults will be saved to path:\n\t\t" << path_output << "/"
              << std::endl;
    std::cout << "\tConfigurations will be read from:\n\t\t" << path_config << "/"
              << std::endl;
    std::cout << "\tConfigurations will be hyp smeared with parameter set "
              << "(alpha1, alpha2, N):\n\t\t" << alpha1 << ", " << alpha2 << ", "
              << iterations << std::endl;
  } catch (std::exception &e) {
    std::cout << e.what() << "\n";
    exit(1);
  }
}

/** Simplifies and cleans GlobalData::read_parameters()
 *
 *  @todo Does this actually make it easier?
 */
void eigenvec_perambulator_input_data_handling(const int number_of_eigen_vec,
                                               const std::string path_eigenvectors,
                                               const std::string name_eigenvectors,
                                               const std::string path_perambulators,
                                               const std::string name_perambulators) {
  try {
    if (number_of_eigen_vec < 1) {
      std::cout << "\ninput file error:\n"
                << "\toption \"number_of_eigen_vec\""
                << " is mandatory and its value must be an integer greater than 0!"
                << "\n\n";
      exit(1);
    } else {
      std::cout << "\tnumber of eigen vectors .................. " << number_of_eigen_vec
                << "\n";
    }
    std::cout << "\tEigenvectors will be read from files:\n\t\t" << path_eigenvectors
              << "/" << name_eigenvectors << "\".eigenvector.t.config\"\n\n";
  } catch (std::exception &e) {
    std::cout << e.what() << "\n";
    exit(1);
  }
}

/** Simplifies and cleans GlobalData::read_parameters()
 *
 *  @todo Does this actually make it easier?
 */
void config_input_data_handling(const int start_config,
                                const int end_config,
                                const int delta_config) {
  try {
    if (start_config < 0) {
      std::cout << "\ninput file error:\n"
                << "\toption \"start config\""
                << " is mandatory and its value must be an integer greater or equal 0!"
                << "\n\n";
      exit(1);
    } else if (end_config < 0 || end_config < start_config) {
      std::cout << "\ninput file error:\n"
                << "\toption \"end_config\""
                << " is mandatory, its value must be an integer greater than 0,"
                << " and it must be larger than start config!"
                << "\n\n";
      exit(1);
    } else if (delta_config < 0) {
      std::cout << "\ninput file error:\n"
                << "\toption \"delta_config\""
                << " is mandatory and its value must be an integer greater than 0!"
                << "\n\n";
      exit(1);
    } else
      std::cout << "\tprocessing configurations " << start_config << " to " << end_config
                << " in steps of " << delta_config << "\n\n";
  } catch (std::exception &e) {
    std::cout << e.what() << "\n";
    exit(1);
  }
}

/** Creates quarks from quark_configs and performs sanity check
 *
 *  @param          quark_configs   Quarks as read from the infile
 *                                  @todo rename to quark_string
 *  @param[in,out]  quarks          Quarks munged into a quark struct.
 *
 *  The work of splitting up the strings and extracting the information is
 *  handled by gdu::make_quark(). The checks are delgated to gdu:quark_check()
 */
void quark_input_data_handling(const std::vector<std::string> quark_configs,
                               std::vector<quark> &quarks) {
  try {
    // Transform each configured quark into a quark via make_quark,
    // inserting each object into the quark vector.
    std::transform(quark_configs.begin(),
                   quark_configs.end(),
                   std::back_inserter(quarks),
                   make_quark);
    // setting id's in quarks
    ssize_t quark_counter = 0;
    for (auto &q : quarks) {
      q.id = quark_counter;
      quark_counter++;
    }
    // checking the contents for correctness
    std::for_each(quarks.begin(), quarks.end(), quark_check);
  } catch (std::exception &e) {
    std::cout << e.what() << "\n";
    exit(1);
  }
}

/** Creates operator_list from operator_list configs
 *
 *  @param          operator_list_configs Operators as read from the infile
 *                                        @todo rename to operator_string
 *  @param[in,out]  operator_list         Operators munged into a Operator
 *                                        struct.
 *
 *  The work of splitting up the strings and extracting the information is
 *  handled by gdu::make_operator_list
 */
void operator_input_data_handling(const std::vector<std::string> operator_strings,
                                  Operator_list &operator_list) {
  try {
    // Transform each configured operator into an Operator_list via
    // make_operator_list()
    for (auto operator_string : operator_strings)
      operator_list.push_back(make_operator_list(operator_string));
    // TODO write a check for correctness of input
  } catch (std::exception &e) {
    std::cout << "operator_input_data_handling: " << e.what() << "\n";
    exit(1);
  }
}

/** Creates correlator_list from correlator_string configs
 *
 *  @param          correlator_string @parblock
 *    Correlators as read from the infile
 *
 *    correlator_list = @em type : @em quark : @em operator : ... : [@em GEVP] :
 *                      [@em P]
 *    where the following abbreviationswhere used
 *    - @em type {C1,C2c,C20,C20V,C3c,C30,C4cD,C4cV,C4cC,C4cB,C40D,C40V,C40C,C40B} :
 *                              Identifier for the Wick diagram to be
 *                              calculated. @see { LapH::Correlators }
 *    - @em quark {"Q%d"} :     Specifies which of the quarks from the infile
 *                              to use
 *    - @em operator {"Op%d"} : Specifies which of the operators from the
 *                              infile to use
 *    - @em GEVP                @todo is that even supported?
 *    - @em P                   @todo is that even supported?
 *
 *    The number of quarks and operators to be specified depends on the diagram
 *    chosen.
 *  @endparblock
 *  @param[in,out]  correlator_list   Correlators munged into a Correlators
 *                                    struct.
 */
void correlator_input_data_handling(const std::vector<std::string> &correlator_strings,
                                    Correlator_list &correlator_list) {
  try {
    // Transform each configured correlator into an Correlator_list via
    // make_correlator()
    for (auto correlator_string : correlator_strings) {
      correlator_list.push_back(make_correlator(correlator_string));
    }
  } catch (std::exception &e) {
    std::cout << "correlator_input_data_handling: " << e.what() << "\n";
    exit(1);
  }
}

}  // end of unnamed namespace

/**
 *  @param[out] gd Global program state.
 *  @param[in] quark_configs Quarks as read from the infile.
 *  @param[in] operator_list_configs Operators as read from the infile.
 *  @param[in] correlator_list_configs Correlators as read from the infile.
 *
 *  @todo Split that into multiple functions for checks, output and munging to
 *  improve readability
 */
void input_handling(GlobalData &gd,
                    std::vector<std::string> const &quark_configs,
                    std::vector<std::string> const &operator_list_configs,
                    std::vector<std::string> const &correlator_list_configs) {
  // Checks and terminal output for lattice, config and paths
  lattice_input_data_handling(gd.path_output,
                              gd.name_lattice,
                              gd.path_config,
                              gd.hyp_parameters.alpha1,
                              gd.hyp_parameters.alpha2,
                              gd.hyp_parameters.iterations,
                              gd.Lt,
                              gd.Lx,
                              gd.Ly,
                              gd.Lz);
  config_input_data_handling(gd.start_config, gd.end_config, gd.delta_config);
  eigenvec_perambulator_input_data_handling(gd.number_of_eigen_vec,
                                            gd.path_eigenvectors,
                                            gd.name_eigenvectors,
                                            gd.path_perambulators,
                                            gd.name_perambulators);

  // Munging of quarks, operators und correlators
  quark_input_data_handling(quark_configs, gd.quarks);
  operator_input_data_handling(operator_list_configs, gd.operator_list);
  correlator_input_data_handling(correlator_list_configs, gd.correlator_list);
}
