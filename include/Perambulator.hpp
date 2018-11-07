/*! @file Perambulator.h
 *  Class decleration of LapH::Perambulator
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 */

#pragma once

#include "global_data_typedefs.hpp"
#include "typedefs.hpp"

#include <Eigen/Dense>

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

/*! Memory allocation and reading routines for perambulators */
class Perambulator {
 private:
  std::vector<Eigen::MatrixXcd> peram;

 public:
  /*! Default constructor */
  Perambulator() : peram(0, Eigen::MatrixXcd(0, 0)){};

  /*! Constructor which creates the correct sizes of the matrices
   *
   *  @param nb_entitites Number of perambulators
   *  @param size_rows    Number of rows for each perambulator
   *  @param size_cols    Number of columns for each perambulator
   */
  Perambulator(const ssize_t nb_entities,
               const std::vector<ssize_t> &size_rows,
               const std::vector<ssize_t> &size_cols)
      : peram(nb_entities, Eigen::MatrixXcd(0, 0)) {
    // TODO: Think about putting this in initialisation list (via lambda?)
    for (ssize_t i = 0; i < nb_entities; i++)
      peram[i].resize(size_rows[i], size_cols[i]);

    std::cout << "\tPerambulators initialised" << std::endl;
  }

  /*! Default deconstructor - std::vector and Eigen should handle everything */
  ~Perambulator(){};

  /*! Overloading [] operator for Perambulator objects */
  inline const Eigen::MatrixXcd &operator[](const ssize_t entity) const {
    return peram.at(entity);
  }

  /*! Reading one perambulators from a single file */
  void read_perambulator(const ssize_t entity,
                         const ssize_t Lt,
                         const ssize_t nb_eigen_vec,
                         const quark &quark,
                         const std::string &filename,
                         const bool mock = false);

  /*! Reading perambulators where each perambulator is stored in a different
   *  file
   */
  void read_perambulators_from_separate_files(
      const ssize_t Lt,
      const ssize_t nb_eigen_vec,
      const std::vector<quark> &quark,
      const std::vector<std::string> &filename_list);
};
