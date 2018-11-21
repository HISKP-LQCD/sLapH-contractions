#pragma once

#include "typedefs.hpp"

#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <string>

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class EigenVector {
 private:
  vec_Xcd_eigen V;

 public:
  // standard ctor that does nothing
  EigenVector() : V(0, Eigen::MatrixXcd(0, 0)){};
  // ctor that creates memory space for the Eigenvectors
  // input: length   -> how many 'timeslices' of eigenvectors are needed
  //        row_size -> length of each eigenvector
  //        col_size -> total number of eigenvectors
  EigenVector(const ssize_t length, const ssize_t row_size, const ssize_t col_size)
      : V(length, Eigen::MatrixXcd(row_size, col_size)){};
  // standard dtor - std::vector and Eigen should be able to handle everything
  ~EigenVector(){};

  // [] operator to directly access the elements of V
  // input: t -> the 'timeslace' from which you want to read the eigenvectors
  inline const Eigen::MatrixXcd &operator[](ssize_t t) const { return V[t]; }

  // writing the eigen vector to some file
  // input: filename -> the filename
  // TODO: NOT IMPLEMENTED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  void write_eigen_vector(const std::string &filename,
                          const ssize_t,
                          const ssize_t verbose);
  // reading the eigen vectors from some file
  // input: filename -> the filename which must contain the path and the
  //                    filename WITHOUT the timeslice (do not forget the dot :)
  //        verbose  -> if 1 additional information will be written out
  void read_eigen_vector(const std::string &filename, const ssize_t verbose);
  // reading the eigen vector from some file for a specific timeslice
  // input: filename -> the path with the FULL filename
  //        t        -> timeslice in V where the eigenvector will be written to
  //        verbose  -> if 1 additional information will be written out
  void read_eigen_vector(const std::string &filename,
                         const ssize_t t,
                         const ssize_t verbose,
                         const bool mock = false);
  void set_V(Eigen::MatrixXcd &v, const ssize_t t);
};
