#pragma once

#include "ComplexProduct.hpp"
#include "OperatorsForMesons.hpp"
#include "global_data.hpp"
#include "timings.hpp"
#include "typedefs.hpp"

#include <omp.h>
#include <Eigen/Dense>

#include <iosfwd>
#include <set>
#include <sstream>
#include <vector>

typedef uint16_t UsedRnd;

struct DilutedFactor {
  using Data = Eigen::MatrixXcd;

  Data data;
  std::pair<RndId, RndId> ric;
  UsedRnd used_rnd_ids;
};

/** Product yielding the off-diagonal elements.

  From the sets of DilutedFactor elements, the product set of DilutedFactor is build such
  that it only contains elements with _unequal_ left and right random vector index. This
  set is intended to be used as an intermediate result.
  */
std::vector<DilutedFactor> operator*(std::vector<DilutedFactor> const &left_vec,
                                     std::vector<DilutedFactor> const &right_vec);

template <int n>
using DilutedFactorsMap = std::map<std::array<ssize_t, n>, std::vector<DilutedFactor>>;

template <size_t n>
std::string to_string(typename DilutedFactorsMap<n>::key_type const &array) {
  std::ostringstream oss;
  oss << "{";
  for (int i = 0; i < n; ++i) {
    if (i != 0) {
      oss << ", ";
    }
    oss << array[i];
  }
  oss << "}";

  return oss.str();
}

template <size_t n>
void print(DilutedFactorsMap<n> const &otfm) {
  std::cout << "DilutedFactorsMap, size = " << otfm.size() << "\n";
  for (auto const &elem : otfm) {
    std::cout << "  " << to_string<n>(elem.first) << " -> "
              << "std::vector(size = " << elem.second.size() << ")\n";
  }
}

Complex trace(std::vector<DilutedFactor> const &left_vec,
              std::vector<DilutedFactor> const &right_vec);

#if 0
int constexpr max_flavor = 8;
using UpTo = boost::container::static_vector<RndId, max_flavor>;

UpTo get_max_used(DilutedTrace const &df, UpTo const &rnd_offset) {
  UpTo result(rnd_offset.size() - 1, 0);

  // Iterate through every random vector index that is either used internally or
  // externally.
  for (auto const id : df.used_rnd_ids) {
    // Iterate through the quark flavors that we have.
    for (int f = 0; f < rnd_offset.size() - 1; ++f) {
      // Does the current random vector index belong to the current flavor?
      if (rnd_offset[f] <= id && id < rnd_offset[f + 1]) {
        // Is this random vector index larger than the largest one stored so far? If so,
        // update this.
        if (result[f] < id) {
          result[f] = id;
        }
      }
    }
  }

  return result;
}

std::map<UpTo, Complex> sub_accumulate(GlobalData const &gd,
                                       std::vector<DilutedTrace> const &traces) {
  // Extract the number of random vectors per quark flavor.
  auto const &quarks = gd.quarks;
  UpTo rnd_count;
  std::transform(std::begin(quarks),
                 std::end(quarks),
                 std::back_inserter(rnd_count),
                 [](quark const &q) { return q.number_of_rnd_vec; });

  // Obtain the offsets for the random vector ids per quark flavor. This basically is an
  // exclusive scan, but this is not in the standard library until C++17.
  UpTo rnd_offset;
  rnd_offset.push_back(0);
  std::partial_sum(
      std::begin(rnd_count), std::end(rnd_count), std::back_inserter(rnd_offset));

  std::map<UpTo, Complex> results;

  // We start with the case where in every flavor we only use 1 random vector.
  UpTo upto(quarks.size(), 1);

  // We want to go until we have used all available random vectors.
  UpTo upto_max = rnd_count;

  // Accumulate the individual contributions into a bucket labeled with the highest random
  // vector index used per flavor.
  for (auto const &trace : traces) {
    auto const &max_used = get_max_used(trace, rnd_offset);
    results[max_used] += trace.data;
  }

  // In the result we want to also accumulate everything below into the current one. There
  // might be an efficient n-dimensional scan algorithm, but this will very likely not be
  // the performance critical part, therefore we just implement it in an easy way with a
  // deep copy.
  auto const parts = results;

  for (auto &result : results) {
    int summands = 0;
    for (auto const &part : parts) {
      // We must skip the element itself, otherwise it would be accumulated onto itself.
      if (result.first == part.first) {
        continue;
      }

      // Take the difference in the maximum random vector index for each flavor.
      auto diff = result.first;
      for (int i = 0; i < diff.size(); ++i) {
        diff[i] -= part.first[i];
      }

      // Are none the differences positive? This means that all random vector indices used
      // are smaller or equal and therefore this contributed. The case where all
      // differences are zero is excluded already at this point.
      bool const none_positive =
          std::none_of(std::begin(diff), std::end(diff), [](RndId d) { return d > 0; });
      // If so, accumulate this onto the result.
      if (none_positive) {
        result.second += part.second;
        ++summands;
      }
    }

    // Normalize.
    result.second /= summands;
  }

  return results;
}
#endif

template <int n1, int n2>
void multiply(DilutedFactorsMap<n1 + n2> &L,
              std::array<ssize_t, n1 + n2> const &key,
              DilutedFactorsMap<n1> const &factor0,
              DilutedFactorsMap<n2> const &factor1) {
  if (L.count(key) == 0) {
    TimingScope<5> timing_scope("multiply");

    // Extract quantum number keys.
    std::array<ssize_t, n1> key1;
    std::array<ssize_t, n2> key2;
    std::copy_n(std::begin(key) + 0, n1, std::begin(key1));
    std::copy_n(std::begin(key) + n1, n2, std::begin(key2));

    auto const &f0 = factor0.at(key1);
    auto const &f1 = factor1.at(key2);

    L[key] = f0 * f1;
  }
}

std::ostream &operator<<(std::ostream &os, DilutedFactor const &df);
