#include "RandomVector.hpp"

#include "typedefs.hpp"

#include "ranlxs.h"

void RandomVector::set(const ssize_t entity, const int seed) {
  // initialisation of the rando vector to create Z2 random vector
  rlxs_init(0, seed);
  std::vector<float> rnd(2 * length);
  ranlxs(&(rnd[0]), 2 * length);

  // generating a Z_2 source
  for (ssize_t i = 0; i < length; ++i) {
    const double sqrt2 = 0.5 * sqrt(2.0);
    double re, im;
    if (rnd[2 * i] < 0.5)
      re = -sqrt2;
    else
      re = sqrt2;
    if (rnd[2 * i + 1] < 0.5)
      im = -sqrt2;
    else
      im = sqrt2;
    vec[entity * length + i] = Complex(re, im);
  }
}

void RandomVector::set(const ssize_t entity,
                       const int seed,
                       const std::string &filename) {
  set(entity, seed);
  write_random_vector(filename);
}

void RandomVector::write_random_vector(const std::string &filename) const {
  // writing random vector to file
  FILE *fp = NULL;
  if ((fp = fopen(filename.c_str(), "wb")) == NULL) {
    std::ostringstream oss;
    oss << "failed to open file to write random vector: " << filename;
    std::runtime_error(oss.str());
  }
  int check_read_in = fwrite(&(vec[0]), sizeof(Complex), vec.size(), fp);
  if (check_read_in != ssize(vec))
    std::cout << "It seems that not all data were written to: " << filename.c_str()
              << "\n"
              << std::endl;
  fclose(fp);
}

void RandomVector::write_random_vector(const ssize_t entity,
                                       const std::string &filename) const {
  // writing random vector to file
  FILE *fp = NULL;
  if ((fp = fopen(filename.c_str(), "wb")) == NULL) {
    std::ostringstream oss;
    oss << "failed to open file to write random vector: " << filename;
    std::runtime_error(oss.str());
  }
  int check_read_in = fwrite(&(vec[entity * length]), sizeof(Complex), length, fp);
  if (check_read_in != length)
    std::cout << "It seems that not all data were written to: " << filename.c_str()
              << "\n"
              << std::endl;
  fclose(fp);
}

void RandomVector::read_random_vector(const std::string &filename) {
  // open file for reading
  FILE *fp = NULL;
  if ((fp = fopen(filename.c_str(), "rb")) == NULL) {
    std::ostringstream oss;
    oss << "failed to open file to read random vector: " << filename;
    std::runtime_error(oss.str());
  }
  // set random vector to zero
  std::fill(vec.begin(), vec.end(), Complex(.0, .0));
  // reading data
  int check_read_in = fread(&(vec[0]), sizeof(Complex), vec.size(), fp);
  if (check_read_in != ssize(vec))
    std::cout << "It seems that not all data are written to: " << filename.c_str() << "\n"
              << std::endl;
  fclose(fp);
}

void RandomVector::read_random_vector(const ssize_t entity, const std::string &filename) {
  // open file for reading
  std::cout << "\tReading random vector from file:\n\t\t" << filename << std::endl;
  FILE *fp = NULL;
  if ((fp = fopen(filename.c_str(), "rb")) == NULL) {
    std::ostringstream oss;
    oss << "failed to open file to read random vector: " << filename;
    std::runtime_error(oss.str());
  }
  // reading data
  int check_read_in = fread(&(vec[entity * length]), sizeof(Complex), length, fp);
  if (check_read_in != length)
    std::cout << "It seems that not all data are read from: " << filename.c_str() << "\n"
              << std::endl;
  fclose(fp);
}

void RandomVector::read_random_vectors_from_separate_files(
    const std::vector<std::string> &filename_list) {
  if (ssize(filename_list) != nb_entities)
    std::cout << "Problem when reading random vectors: The number of random "
              << "vectors read is not the same as the expected one!" << std::endl;
  // set random vector to zero
  std::fill(vec.begin(), vec.end(), Complex(.0, .0));
  for (ssize_t i = 0; i < ssize(filename_list); i++)
    read_random_vector(i, filename_list[i]);
}
