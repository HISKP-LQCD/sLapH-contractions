#pragma once

#include "typedefs.hpp"

#include <ostream>

/*! Small struct which contains all information to build and read random
 *  vectors.
 *
 *  How many are needed, their length, and what the input paths are.
 */
struct RandomVectorConstruction {
  ssize_t nb_entities, length;
  std::vector<std::string> filename_list;
};

/*! Small struct which contains all information to build and read
 *  Perambulators.
 *
 *  How many are needed, their size, and what the input paths are.
 */
struct PerambulatorConstruction {
  ssize_t nb_entities;
  std::vector<ssize_t> size_rows;
  std::vector<ssize_t> size_cols;
  std::vector<std::string> filename_list;
};

/*! Quark type that contains all quark propagator informations:
 *
 *  Flavor, number of random vectors, dilution scheme, path and a unique id.
 */
struct quark {
  std::string type; /*!< Flavor */
  int number_of_rnd_vec;

  std::string dilution_T;
  int number_of_dilution_T;

  std::string dilution_E;
  int number_of_dilution_E;

  std::string dilution_D;
  int number_of_dilution_D;

  ssize_t id;
  std::string path;

  /*! Constructor */
  quark(std::string type,
        int number_of_rnd_vec,
        std::string dilution_T,
        int number_of_dilution_T,
        std::string dilution_E,
        int number_of_dilution_E,
        std::string dilution_D,
        int number_of_dilution_D,
        int id,
        std::string path)
      : type(type),
        number_of_rnd_vec(number_of_rnd_vec),
        dilution_T(dilution_T),
        number_of_dilution_T(number_of_dilution_T),
        dilution_E(dilution_E),
        number_of_dilution_E(number_of_dilution_E),
        dilution_D(dilution_D),
        number_of_dilution_D(number_of_dilution_D),
        id(id),
        path(path) {}
};

/*! Struct that contains all physical information specifying a quantum field
 *  operator in the infile:
 *
 *  Dirac structure, displacement and momentum vector as 3-vectors
 *
 *  @todo rewrite the momenta with eigen or at least overload +, - and abs for
 *        them
 */
struct QuantumNumbers {
  using VectorData = Eigen::Vector3i;
  //  using VectorData = std::array<int, 3>;

  std::vector<int> gamma;
  DisplacementDirection displacement;
  VectorData momentum;
};

inline std::string to_string(std::vector<int> const &vec) {
  std::string result("");
  for (int el : vec) {
    result += std::to_string(el);
  }
  return result;
}

inline std::string to_string(Eigen::Vector3i const &vec) {
  std::string result("");
  for (int i = 0; i < 3; i++) {
    result += std::to_string(vec[i]);
  }
  return result;
}

/*! Makes a string object of a displacement vector */
inline std::string to_string(DisplacementDirection const &vec){
  std::string out;
  if (vec.empty()) out = "000";
  for (auto const& dis : vec){ 
    out.push_back(dis.first);
    out.push_back(dis.second);
  }
  return out;
}
inline std::ostream &operator<<(std::ostream &os, QuantumNumbers const &qn) {
  os << "\tmomentum: " << qn.momentum[0] << qn.momentum[1] << qn.momentum[2]
     << "\n\tdisplacement: ";
  for (auto const &dir : qn.displacement) {
    os << dir.first << dir.second;
  }
  os << "\n\tgamma struct: ";
  for (const auto &g : qn.gamma)
    os << g;
  os << "\n" << std::endl;

  return os;
}

/*! Struct that contains all information specifying the correlator in the
 *  infile:
 *
 *  Internal name of the diagram (@ref LapH::Correlators), and which quark from
 *  quarks and operator from oeprator_list are to be used
 *
 *  @todo Are GEVP and tot_mom used at the moment? They should support
 *        building off-diagonal elements as well or resticting to a certain
 *        CMS momentum
 */
struct Correlators_2 {
  using VectorData = Eigen::Vector3i;

  std::string type;
  std::vector<int> quark_numbers;
  std::vector<int> operator_numbers;
  std::string GEVP;
  std::vector<VectorData> tot_mom;
};

typedef std::vector<QuantumNumbers> Operators;
typedef std::vector<Operators> Operator_list;
typedef std::vector<Correlators_2> Correlator_list;

//Typedef for hypercubic blocking of gauge links, takes alpha1, alpha2 and
//number of iterations
struct HypPars {
  double alpha1;
  double alpha2;
  ssize_t iterations;
};
