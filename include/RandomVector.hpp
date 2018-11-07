#pragma once

#include "typedefs.hpp"

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/** Memory allocation and IO routines for random vectors */
class RandomVector {
 private:
  // the random vector
  std::vector<Complex> vec;
  const ssize_t nb_entities;
  const ssize_t length;

 public:
  /** Constructor which delegates memory allocation to std::vector
   *
   *  @param nb_entities  Number of random vectors
   *  @param length       Number of complex entries in each random vector
   *
   *  @warning It is assumed that all random vector have the same length!
   */
  RandomVector(const ssize_t nb_entities, const ssize_t length)
      : vec(nb_entities * length, Complex(0.0, 0.0)),
        nb_entities(nb_entities),
        length(length) {
    std::cout << "\tRandom vectors initialised" << std::endl;
  };

  /** Default deconstructor */
  ~RandomVector(){};

  /** Overloading () operator to directly access the elements of vec */
  inline Complex operator()(const ssize_t entity, const ssize_t entry) const {
    return vec.at(entity * length + entry);
  }

  const std::vector<Complex> & get(void) const {
    return vec;
  }

  /** Computes the random vectors for the sources
   *
   *  @param entity The random vector you want to set
   *  @param seed   Seed for the random vector
   *  @param length Length of the random vector
   */
  void set(const ssize_t entity, const int seed);

  /** Computes the random vectors for the sources and stores them to file
   *
   *  @param entity   The random vector you want to set
   *  @param seed     Seed for the random vector
   *  @param length   Length of the random vector
   *  @param filename The random vector is directly stored in this file
   */
  void set(const ssize_t entity, const int seed, const std::string &filename);

  /** Write all random vectors to some file */
  void write_random_vector(const std::string &filename) const;

  /** Write one random vector to some file
   *
   *  @param entity   The random vector to be stored on disk
   *  @param filename The filename
   */
  void write_random_vector(const ssize_t entity, const std::string &filename) const;

  /** Read all random vectors from some file */
  void read_random_vector(const std::string &filename);

  /** Read one random vector from some file
   *
   *  @param entity   The random vector to be read from disk
   *  @param filename The filename
   */
  void read_random_vector(const ssize_t entity, const std::string &filename);

  /** Read random vectors where each vector is stored in a different file
   *
   *  @param filename_list Vector which contains all file names
   */
  void read_random_vectors_from_separate_files(
      const std::vector<std::string> &filename_list);
};
