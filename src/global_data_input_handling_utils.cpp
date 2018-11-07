#include "global_data.hpp"
#include "global_data_utils.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>


namespace {

using Vector = QuantumNumbers::VectorData;

// *****************************************************************************
// A helper function to simplify the main part.
template <class T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
  return os;
}
// *****************************************************************************
/// @brief Stream insertion operator for slave.
///
/// @param stream The stream into which quark is being inserted.
/// @param q The quark object.
///
/// @return Reference to the ostream.
inline std::ostream &operator<<(std::ostream &stream, const quark &quark) {
  return stream << "\tQUARK type: ****  " << quark.type
                << "  ****\n\t number of random vectors: " << quark.number_of_rnd_vec
                << "\n\t dilution scheme in time: " << quark.dilution_T
                << quark.number_of_dilution_T
                << "\n\t dilution scheme in ev space: " << quark.dilution_E
                << quark.number_of_dilution_E
                << "\n\t dilution scheme in Dirac space: " << quark.dilution_D
                << quark.number_of_dilution_D
                << "\n\t path of the perambulator and random vectors:\n\t\t" << quark.path
                << "\n";
}

// *****************************************************************************
inline Vector create_3darray_from_string(std::string in) {
  Vector out;
  std::vector<std::string> tokens;
  // erasing the brakets at the beginning and the end
  in.erase(0, 2);
  in.erase(in.end() - 1);

  boost::split(tokens, in, boost::is_any_of(","));

  return {boost::lexical_cast<int>(tokens[0]),
          boost::lexical_cast<int>(tokens[1]),
          boost::lexical_cast<int>(tokens[2])};
}
// *****************************************************************************
inline void create_all_momentum_combinations(const int p, std::vector<Vector> &out) {
  // creating all momentum combinations possible and needed
  int max_p = p;
  std::vector<Vector> all_p;
  for (int p1 = -max_p; p1 < max_p + 1; p1++)
    for (int p2 = -max_p; p2 < max_p + 1; p2++)
      for (int p3 = -max_p; p3 < max_p + 1; p3++)
        all_p.push_back({p1, p2, p3});
  // copying wanted combinations into out array
  for (const auto &all : all_p)
    if (p == all[0] * all[0] + all[1] * all[1] + all[2] * all[2])
      out.push_back(all);
}
// *****************************************************************************
inline void create_mom_array_from_string(std::string in,
                                         std::vector<std::vector<Vector>> &out) {
  // erase the p (first entry)
  in.erase(0, 1);
  std::vector<std::string> tokens;
  boost::split(tokens, in, boost::is_any_of(","));
  int p;
  ssize_t counter = 0;
  out.resize(tokens.size());
  for (const auto &t : tokens) {
    p = boost::lexical_cast<int>(t);
    create_all_momentum_combinations(p, out[counter]);
    counter++;
  }
}

DisplacementDirection make_displacement_direction(std::vector<std::string> token) {
  DisplacementDirection result;
  for (auto &direction : token) {
    result.push_back(std::make_pair(direction[0], direction[1]));
  }

  return result;
}

// *****************************************************************************
inline void create_displacement_direction_from_string(
    std::string in, std::vector<DisplacementDirection> &out) {
  // erase the d (first entry)
  in.erase(0, 1);
  std::vector<std::string> tokens;
  boost::split(tokens, in, boost::is_any_of(","));
  out.resize(tokens.size());
  for (ssize_t counter = 0; counter < ssize(tokens); ++counter) {
    std::vector<std::string> t;
    boost::split(t, tokens[counter], boost::is_any_of("|"));
    out[counter] = make_displacement_direction(t);
  }
}

}  // end of unnamed namespace
/******************************************************************************/
/******************************************************************************/

namespace global_data_utils {

/*!
 *  @param quark_string Quark as specified in the infile
 *                      quark = @em flavor : @em nb_rnd_vec : T @em diltT_type :
 *                              @em dilT : E @em dilE_type : @em dilE :
 *                              D @em dilD_type : @em dilD : @em path
 *  where the following abbreviations where used
 *  - @em flavor     {u,s,c,b} : Quark flavor
 *  - @em nb_rnd_vec       : Number of random vectors
 *  - @em dilT_type  {I,B} : Dilution type in time
 *  - @em dilT             : Number of dilution blocks in time
 *  - @em dilE_type  {I,B} : Dilution type in eigenvetor space
 *  - @em dilE             : Number of dilution blocks in eigenvector space
 *  - @em dilD_type  {I,B} : Dilution type in dirac space
 *  - @em dilD             : Number of dilution blocks in Dirac space
 *  - @em path             : Path to perambulator
 *  The validity of the values is checked in quark_check()
 *
 *  @returns A quark object constructed with the data obtained from
 *           @em quark_string
 *
 *  Internally uses boost to split the string and process the parts.
 */
quark make_quark(const std::string &quark_string) {
  // Tokenize the string on the ":" delimiter.
  std::vector<std::string> tokens;
  boost::split(tokens, quark_string, boost::is_any_of(":"));

  // If the split did not result in exactly 8 tokens, then the value
  // is formatted wrong.
  if (9 != tokens.size()) {
    using boost::program_options::validation_error;
    throw validation_error(
        validation_error::invalid_option_value, "quarks.quark", quark_string);
  }

  // Create a quark from the token values.
  return quark(tokens[0],
               boost::lexical_cast<int>(tokens[1]),
               tokens[2],
               boost::lexical_cast<int>(tokens[3]),
               tokens[4],
               boost::lexical_cast<int>(tokens[5]),
               tokens[6],
               boost::lexical_cast<int>(tokens[7]),
               0,
               tokens[8]);
}

void quark_check(quark quarks) {
  if (quarks.type != "u" && quarks.type != "d" && quarks.type != "s" &&
      quarks.type != "c") {
    throw std::runtime_error("quarks.quark.type must be u, d, s or c");
  }
  /*! @todo Check that the number of random vectors is greater than the
   *  largest required diagram
   */
  else if (quarks.number_of_rnd_vec < 1) {
    throw std::runtime_error("quarks.quark.number_of_rnd_vec must be greater than 0");
  } else if (quarks.dilution_T != "TI" && quarks.dilution_T != "TB" &&
             quarks.dilution_T != "TF") {
    throw std::runtime_error("quarks.quark.dilution_T must be TI, TB, TF");
  } else if (quarks.number_of_dilution_T < 1) {
    throw std::runtime_error(
        "quarks.quark.number_of_dilution_T must be greater than 0 and smaller than the "
        "temporal extent");
  } else if (quarks.dilution_E != "EI" && quarks.dilution_E != "EB" &&
             quarks.dilution_E != "EF") {
    throw std::runtime_error("quarks.quark.dilution_E must be EI, EB or EF");
  } else if (quarks.number_of_dilution_E < 1) {
    throw std::runtime_error(
        "quarks.quark.number_of_dilution_E must be greater than 0 and smaller than "
        "number of eigen vectors");
  } else if (quarks.dilution_D != "DI" && quarks.dilution_D != "DB" &&
             quarks.dilution_D != "DF") {
    throw std::runtime_error("quarks.quark.dilution_D must be DI, DB or DF");
  } else if (quarks.number_of_dilution_D < 1 || quarks.number_of_dilution_D > 4) {
    throw std::runtime_error(
        "quarks.quark.number_of_dilution_D must be greater than 0 and smaller than 5");
  } else {
    std::cout << quarks << std::endl;
  }
}

/*****************************************************************************/
/*!
 *  @param operator_string  Operator as specified in the infile:
 *                          A ';'-sperated list with individual operators. The
 *                          individual operators are composed of '.'-seperated
 *                          parts. E.g.
 *                          @code
 *                            operator_list = g4.d0.p(0,0,1);g5.d0.p0,1
 *                          @endcode
 *                          Momenta (<em>p</em>) can be specified as 3-momentum or
 *                          by (one or more) scalar number(s). In the latter
 *                          case all 3-momenta with corresponding absolute
 *                          value are constructed
 *
 *  @returns An Operator_list object constructed with the data obtained from
 *           @em operator_string
 *
 *  Internally uses boost to split the string and process the parts.
 *
 *  This is a factory function which returns by value and uses named return
 *  value optimization in order to call the constructor within the scope
 */
Operators make_operator_list(const std::string &operator_string) {
  Operators op_list;  // return object

  // Two steps are necessary:
  // 1. Getting all operators in one list which are separated by ";"
  // 2. Separating the individual operators into its smaller bits, which are
  //    separated by "."
  // Tokenize the string on the ";" delimiter -> Individual operators
  std::vector<std::string> operator_tokens;
  boost::split(operator_tokens, operator_string, boost::is_any_of(":"));

  // running over opeator tokens and split them further (Step 2):
  for (const auto &op_t : operator_tokens) {
    std::vector<std::string> tokens;
    boost::split(tokens, op_t, boost::is_any_of("."));
    std::vector<int> gammas;
    std::vector<DisplacementDirection> disp_dirs;
    std::vector<std::vector<Vector>> mom_vec;
    for (auto &str : tokens) {
      // getting the gamma structure
      if (str.compare(0, 1, "g") == 0)
        gammas.push_back(boost::lexical_cast<int>(str.erase(0, 1)));
      // getting the displacement indices
      else if (str.compare(0, 1, "d") == 0) {
        // 0 encodes no displacement. disp_dir remains empty in this case
        if (str.compare(1, 1, "0") == 0) {
          disp_dirs.push_back(DisplacementDirection());
        } else if ((str.compare(1, 1, "<") == 0) || (str.compare(1, 1, ">") == 0))
          create_displacement_direction_from_string(str, disp_dirs);
        else {
          throw std::runtime_error("Something wrong with the displacement in the operator definition!");
        }
      }
      // getting the momenta
      else if (str.compare(0, 1, "p") == 0) {
        if (str.compare(1, 1, "(") == 0) {
          mom_vec.resize(1);
          mom_vec[0].push_back(create_3darray_from_string(str));
        } else {
          create_mom_array_from_string(str, mom_vec);
        }
      }
      // catching wrong entries
      else {
        throw std::runtime_error("There is something wrong with the operators!");
      }
    }

    for (const auto &mom_vec_tmp : mom_vec) {  // momenta
      for (auto mom : mom_vec_tmp) {
        for (auto &disp_dir : disp_dirs) {
          op_list.push_back({gammas, disp_dir, mom});
        }
      }
    }
  }

  std::cout << "op_list has size: " << op_list.size() << std::endl;
  return op_list;
}

/*!
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
 *
 *  Internally uses boost to split the string and process the parts.
 *
 *  @todo Write check for correctness of correlator_string
 */
Correlators_2 make_correlator(const std::string &correlator_string) {
  std::vector<std::string> correlator_tokens;
  boost::split(correlator_tokens, correlator_string, boost::is_any_of(":"));

  std::string type;
  std::vector<int> quark_number;
  std::vector<int> operator_number;
  std::string GEVP;
  std::vector<Vector> tot_mom;

  for (auto corr_t : correlator_tokens) {
    // getting the type name
    if (corr_t.compare(0, 1, "C") == 0)
      type = corr_t;
    // getting quark numbers
    else if (corr_t.compare(0, 1, "Q") == 0)
      quark_number.push_back(boost::lexical_cast<int>(corr_t.erase(0, 1)));
    // getting operator numbers
    else if (corr_t.compare(0, 2, "Op") == 0)
      operator_number.push_back(boost::lexical_cast<int>(corr_t.erase(0, 2)));
    // getting the GEVP type
    else if (corr_t.compare(0, 1, "G") == 0)
      GEVP = corr_t;
    // getting total momenta for moving frames
    else if (corr_t.compare(0, 1, "P") == 0) {
      if (corr_t.compare(1, 1, "(") == 0) {
        tot_mom.push_back(create_3darray_from_string(corr_t));
      } else {
        corr_t.erase(0, 1);
        std::vector<std::string> tokens;
        boost::split(tokens, corr_t, boost::is_any_of(","));
        for (auto token : tokens) {
          create_all_momentum_combinations(boost::lexical_cast<int>(token), tot_mom);
        }
      }
    }

    // catching wrong entries
    else {
      throw std::runtime_error("There is something wrong with the correlators in the input file!");
    }
  }

  return {type, quark_number, operator_number, GEVP, tot_mom};
}

}  // end of namespace global_data_utils
