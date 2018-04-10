#pragma once

#include "DilutedFactor.h"
#include "Gamma.h"
#include "OperatorsForMesons.h"
#include "Perambulator.h"
#include "dilution-iterator.h"
#include "typedefs.h"

#include "Eigen/Dense"
#include "boost/circular_buffer.hpp"
#include "boost/multi_array.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

template <typename T, ssize_t n>
void print(std::array<T, n> const &a) {
  for (auto const &elem : a) {
    std::cout << elem << "\t";
  }
  std::cout << std::endl;
}

template <DilutedFactorType qlt>
class DilutedFactorFactory {
 public:
  using Key = std::array<int, DilutedFactorTypeTraits<qlt>::num_times>;
  using Value = DilutedFactors<1, 0>;

  DilutedFactorFactory(
      RandomVector const &random_vector,
      Perambulator const &perambulator,
      OperatorFactory const &_meson_operator,
      ssize_t const dilT,
      ssize_t const dilE,
      ssize_t const nev,
      typename DilutedFactorTypeTraits<qlt>::type const &quarkline_indices);

  Value const &operator[](Key const &key) {
    if (Ql.count(key) == 0) {
      build(key);
    }

    return Ql.at(key);
  }

  void clear() { Ql.clear(); }

  void build(Key const &time_key);

 private:
  std::map<Key, Value> Ql;

  Perambulator const &peram;
  RandomVector const &rnd_vec;
  OperatorFactory const &meson_operator;
  const ssize_t dilT, dilE, nev;
  typename DilutedFactorTypeTraits<qlt>::type const &quarkline_indices;

  static int constexpr dilD = 4;
};
