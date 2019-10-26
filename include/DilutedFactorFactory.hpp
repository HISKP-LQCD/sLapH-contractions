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

  void request(Key const &time_key) { requests_.insert(time_key); }

  Value const &operator[](Key const &time_key) {
    return Ql.at(time_key);
  }

  void build_all() {
    // The requests have been stored in a set. We need to convert them into a vector such
    // that we can run them concurrently.
    std::vector<Key> unique_requests;
    unique_requests.reserve(requests_.size());
    for (auto const &time_key : requests_) {
      // The map might already contain some elements that have been requested in an
      // earlier iteration. Therefore we need to see whether it has already been built
      // before.
      if (Ql.count(time_key) == 0) {
        unique_requests.push_back(time_key);

        // Populate the whole map with all the keys that are going to be built next. This
        // way the map does not change any more and concurrent read access is possible.
        Ql[time_key];
      }
    }

    // We are going to build all the ones that are requested, therefore the list of
    // requests has to be cleared.
    requests_.clear();

    // Build all the elements. The `build` function will automatically populate the map
    // `Ql`.
    for (auto i = 0; i < ssize(unique_requests); ++i) {
      build(unique_requests[i]);
    }
  }

  void clear() { Ql.clear(); }

 private:
  void build(Key const &time_key);

  std::map<Key, Value> Ql;
  std::set<Key> requests_;

  Perambulator const &peram;
  RandomVector const &rnd_vec;
  OperatorFactory const &meson_operator;
  const ssize_t dilT, dilE, nev;
  typename DilutedFactorTypeTraits<qlt>::type const &quarkline_indices;

  static int constexpr dilD = 4;
};
