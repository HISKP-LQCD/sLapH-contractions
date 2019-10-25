#pragma once

#include "DilutedFactor.hpp"
#include "Gamma.hpp"
#include "OperatorsForMesons.hpp"
#include "Perambulator.hpp"
#include "dilution-iterator.hpp"
#include "typedefs.hpp"

#include <Eigen/Dense>
#include <boost/circular_buffer.hpp>
#include <boost/multi_array.hpp>

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
  using Value = DilutedFactorsMap<1>;

  DilutedFactorFactory(
      RandomVector const &random_vector,
      Perambulator const &perambulator,
      OperatorFactory const &_meson_operator,
      ssize_t const dilT,
      ssize_t const dilE,
      ssize_t const nev,
      typename DilutedFactorTypeTraits<qlt>::type const &quarkline_indices);

  Value const &operator[](Key const &time_key) { return Ql.at(time_key); }

  void build_all() {
    // Populate the whole map with all the keys that are going to be built next. This way
    // the map does not change any more and concurrent read access is possible.
    for (auto const time_key : requests_) {
      Ql[time_key];
    }

    // Build all the elements. The `build` function will automatically populate the map
    // `Ql`.
    for (auto i = 0; i < ssize(requests_); ++i) {
      auto const &time_key = requests_[i];
      build(time_key);
    }
  }

  void clear() { Ql.clear(); }

  void request(Key const &time_key) { requests_.push_back(time_key); }

 private:
  void build(Key const &time_key);

  std::map<Key, Value> Ql;
  std::vector<Key> requests_;

  Perambulator const &peram;
  RandomVector const &rnd_vec;
  OperatorFactory const &meson_operator;
  const ssize_t dilT, dilE, nev;
  typename DilutedFactorTypeTraits<qlt>::type const &quarkline_indices;

  static int constexpr dilD = 4;
};
