#pragma once

#include "OperatorsForMesons.hpp"
#include "global_data.hpp"
#include "local_timer.hpp"
#include "typedefs.hpp"

#include <omp.h>
#include <Eigen/Dense>

#include <iosfwd>
#include <set>
#include <sstream>
#include <vector>

template <size_t rvecs1, size_t rvecs2>
bool has_intersection(SmallVectorRndId<rvecs1> const &left,
                      SmallVectorRndId<rvecs2> const &right) {
  SmallVectorRndId<rvecs1 + rvecs2> intersection;

  std::set_intersection(std::begin(left),
                        std::end(left),
                        std::begin(right),
                        std::end(right),
                        std::back_inserter(intersection));
  return intersection.size() > 0;
}

template <size_t rvecs1, size_t rvecs2>
void merge_append(SmallVectorRndId<rvecs1> &data,
                  SmallVectorRndId<rvecs2> const &addition) {
  auto const old_end = std::end(data);
  std::copy(std::begin(addition), std::end(addition), std::back_inserter(data));
  std::inplace_merge(std::begin(data), old_end, std::end(data));
}

template <size_t rvecs>
void merge_push_back(SmallVectorRndId<rvecs> &data, RndId const &addition) {
  auto const old_end = std::end(data);
  data.push_back(addition);
  std::inplace_merge(std::begin(data), old_end, std::end(data));
}

/** Product yielding the off-diagonal elements.

  From the sets of DilutedFactor elements, the product set of DilutedFactor is build such
  that it only contains elements with _unequal_ left and right random vector index. This
  set is intended to be used as an intermediate result.
  */
template <size_t rvecs1, size_t rvecs2>
std::vector<DilutedFactor<rvecs1 + rvecs2 + 1>> operator*(
    std::vector<DilutedFactor<rvecs1>> const &left_vec,
    std::vector<DilutedFactor<rvecs2>> const &right_vec) {
  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  int constexpr rvecs_total = rvecs1 + rvecs2 + 1;
  std::vector<DilutedFactor<rvecs_total>> result_vec;

  for (auto const &left : left_vec) {
    auto const inner_rnd_id = left.ric.second;

    for (auto const &right : right_vec) {
      // We want to make the inner and outer indices differ. The inner indices need to
      // match because the product would not make sense otherwise.
      bool const is_allowed =
          inner_rnd_id == right.ric.first && left.ric.first != right.ric.second;
      if (!is_allowed) {
        continue;
      }

      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if (has_intersection(left.used_rnd_ids, right.used_rnd_ids)) {
        continue;
      }

      // We want to keep track of the indices that have been contracted away. These are
      // all the ones from the left factor, all the ones from the right factor and the one
      // that we are contracting over right now.
      SmallVectorRndId<rvecs_total> used;
      used.push_back(inner_rnd_id);
      merge_append(used, left.used_rnd_ids);
      merge_append(used, right.used_rnd_ids);

      result_vec.push_back({left.data * right.data,
                            std::make_pair(left.ric.first, right.ric.second),
                            used});
    }
  }

  if (result_vec.size() == 0) {
    throw std::runtime_error(
        "vector<DilutedFactor> operator*(vector<DilutedFactor>, vector<DilutedFactor>) "
        "has an empty result.");
  }

  return result_vec;
}

template <size_t rvecs>
inline Complex operator+(DilutedTrace<rvecs> const &df, Complex const &c) {
  return c + df.data;
}

template <size_t rvecs>
inline Complex operator+(Complex const &c, DilutedTrace<rvecs> const &df) {
  return df + c;
}

template <int n, size_t rvecs>
using DilutedFactors =
    std::map<std::array<ssize_t, n>, std::vector<DilutedFactor<rvecs>>>;

template <size_t n, size_t rvecs>
std::string to_string(typename DilutedFactors<n, rvecs>::key_type const &array) {
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

template <size_t n, size_t rvecs>
void print(DilutedFactors<n, rvecs> const &otfm) {
  std::cout << "DilutedFactors, size = " << otfm.size() << "\n";
  for (auto const &elem : otfm) {
    std::cout << "  " << to_string<n, rvecs>(elem.first) << " -> "
              << "std::vector(size = " << elem.second.size() << ")\n";
  }
}

template <size_t rvecs1, size_t rvecs2>
Complex trace(std::vector<DilutedFactor<rvecs1>> const &left_vec,
              std::vector<DilutedFactor<rvecs2>> const &right_vec) {
  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  Complex result(0.0, 0.0);

  int num_summands = 0;

  LT_ULTRA_FINE_DECLARE;

  for (auto const &left : left_vec) {
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

    LT_ULTRA_FINE_START;
    Eigen::MatrixXcd right_sum(
        Eigen::MatrixXcd::Zero(left.data.rows(), left.data.cols()));

    for (auto const &right : right_vec) {
      // We want to make the inner and outer indices match. The inner indices need to
      // match because the product would not make sense otherwise. The outer indices must
      // match since we want to be able to take the trace over the result. The second
      // condition is where this differs from the other multiplication operator.
      bool const is_allowed =
          inner_rnd_id == right.ric.first && outer_rnd_id == right.ric.second;
      if (!is_allowed) {
        continue;
      }

      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if (has_intersection(left.used_rnd_ids, right.used_rnd_ids)) {
        continue;
      }

      // The right sides that we encounter at this point have the same left and right
      // random vector indices. They may differ in the set of used random vector indices.
      // But since we do not plan to contract the result with more DilutedFactor
      // instances, we do not care to preserve the information about the used random
      // vector indices. Therefore we can sum all these elements up to have less
      // multiplications to do.
      right_sum += right.data;
      ++num_summands;
    }
    LT_ULTRA_FINE_STOP;
    LT_ULTRA_FINE_PRINT("[DilutedFactor::trace] right_sum");

    LT_ULTRA_FINE_START;
    auto const &product = left.data * right_sum;
    result += product.trace();
    LT_ULTRA_FINE_STOP;
    LT_ULTRA_FINE_PRINT("[DilutedFactor::trace] product_trace");

  }  // for(left_vec)

  if (num_summands == 0) {
    throw std::runtime_error(
        "cmplx trace(vector<DilutedFactor>, vector<DilutedFactor>) has used zero "
        "summands.");
  }

  return result / static_cast<double>(num_summands);
}

template <size_t rvecs1, size_t rvecs2>
std::vector<DilutedTrace<rvecs1 + rvecs2 + 2>> factor_to_trace(
    std::vector<DilutedFactor<rvecs1>> const &left_vec,
    std::vector<DilutedFactor<rvecs2>> const &right_vec) {
  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  int constexpr rvecs_result = rvecs1 + rvecs2 + 2;
  std::vector<DilutedTrace<rvecs_result>> result_vec;

  LT_ULTRA_FINE_DECLARE;

  for (auto const &left : left_vec) {
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

    LT_ULTRA_FINE_START;
    for (auto const &right : right_vec) {
      // We want to make the inner and outer indices match. The inner indices need to
      // match because the product would not make sense otherwise. The outer indices must
      // match since we want to be able to take the trace over the result. The second
      // condition is where this differs from the other multiplication operator.
      bool const is_allowed =
          inner_rnd_id == right.ric.first && outer_rnd_id == right.ric.second;
      if (!is_allowed) {
        continue;
      }

      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if (has_intersection(left.used_rnd_ids, right.used_rnd_ids)) {
        continue;
      }

      // We want to keep track of the indices that have been contracted away. These are
      // all the ones from the left factor, all the ones from the right factor and the one
      // that we are contracting over right now.
      SmallVectorRndId<rvecs_result> used;
      merge_push_back(used, inner_rnd_id);
      merge_push_back(used, outer_rnd_id);
      merge_append(used, left.used_rnd_ids);
      merge_append(used, right.used_rnd_ids);

      // The right sides that we encounter at this point have the same left and right
      // random vector indices. They may differ in the set of used random vector indices.
      // But since we do not plan to contract the result with more DilutedFactor
      // instances, we do not care to preserve the information about the used random
      // vector indices. Therefore we can sum all these elements up to have less
      // multiplications to do.
      result_vec.push_back(
          {typename DilutedTrace<rvecs_result>::Data{(left.data * right.data).trace()},
           used});
    }
    LT_ULTRA_FINE_STOP;
    LT_ULTRA_FINE_PRINT("[DilutedFactor::factor_to_trace] multiply_trace");
  }  // for(left_vec)

  if (result_vec.size() == 0) {
    throw std::runtime_error(
        "vector<DilutedTrace> factor_to_trace(vector<DilutedFactor>, "
        "vector<DilutedFactor>) has an empty result.");
  }

  return result_vec;
}

template <size_t rvecs>
std::vector<DilutedTrace<rvecs + 1>> factor_to_trace(
    std::vector<DilutedFactor<rvecs>> const &vec) {
  assert(vec.size() > 0);

  std::vector<DilutedTrace<rvecs + 1>> result_vec;

  for (auto const &elem : vec) {
    // We only want to use diagonal elements.
    if (elem.ric.first != elem.ric.second) {
      continue;
    }

    SmallVectorRndId<rvecs + 1> used;
    std::copy(std::begin(elem.used_rnd_ids),
              std::end(elem.used_rnd_ids),
              std::back_inserter(used));

    auto const outer_rnd_id = elem.ric.first;
    merge_push_back(used, outer_rnd_id);

    DilutedTrace<rvecs + 1> result = {elem.data.trace(), used};

    result_vec.push_back(result);
  }

  if (result_vec.size() == 0) {
    throw std::runtime_error(
        "vector<DilutedTrace> factor_to_trace(vector<DilutedFactor>) has an empty "
        "result.");
  }

  return result_vec;
}

template <size_t rvecs1, size_t rvecs2>
ComplexProduct inner_product(std::vector<DilutedTrace<rvecs1>> const &left_vec,
                             std::vector<DilutedTrace<rvecs2>> const &right_vec) {
  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  int num_summands = 0;

  ComplexProduct result(0.0, 0.0, 0.0, 0.0);

  for (auto const &left : left_vec) {
    Complex right_sum(0.0, 0.0);

    for (auto const &right : right_vec) {
      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if (has_intersection(left.used_rnd_ids, right.used_rnd_ids)) {
        continue;
      }

      // The right sides that we encounter at this point have the same left and right
      // random vector indices. They may differ in the set of used random vector indices.
      // But since we do not plan to contract the result with more DilutedFactor
      // instances, we do not care to preserve the information about the used random
      // vector indices. Therefore we can sum all these elements up to have less
      // multiplications to do.
      right_sum += right.data;
      ++num_summands;
    }

    result.rere += left.data.real() * right_sum.real();
    result.reim += left.data.real() * right_sum.imag();
    result.imre += left.data.imag() * right_sum.real();
    result.imim += left.data.imag() * right_sum.imag();
  }

  if (num_summands == 0) {
    throw std::runtime_error(
        "compcomp_t inner_product(vector<DilutedTrace>, "
        "vector<DilutedTrace>) has an empty result.");
  }

  return result / static_cast<double>(num_summands);
}

int constexpr max_flavor = 8;
using UpTo = boost::container::static_vector<RndId, max_flavor>;

template <size_t rvecs>
UpTo get_max_used(DilutedTrace<rvecs> const &df, UpTo const &rnd_offset) {
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

template <size_t rvecs>
std::map<UpTo, Complex> sub_accumulate(GlobalData const &gd,
                                       std::vector<DilutedTrace<rvecs>> const &traces) {
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

template <int n1, int n2, size_t rvecs1, size_t rvecs2>
void multiply(DilutedFactors<n1 + n2, rvecs1 + rvecs2 + 1> &L,
              std::array<ssize_t, n1 + n2> const &key,
              DilutedFactors<n1, rvecs1> const &factor0,
              DilutedFactors<n2, rvecs2> const &factor1) {
  LT_ULTRA_FINE_DECLARE;
  if (L.count(key) == 0) {
    std::array<ssize_t, n1> key1;
    std::array<ssize_t, n2> key2;

    std::copy_n(std::begin(key) + 0, n1, std::begin(key1));
    std::copy_n(std::begin(key) + n1, n2, std::begin(key2));

    auto const &f0 = factor0.at(key1);
    auto const &f1 = factor1.at(key2);

    LT_ULTRA_FINE_START;

    L[key] = f0 * f1;

    LT_ULTRA_FINE_STOP;
    LT_ULTRA_FINE_PRINT("[DilutedFactor::multiply] multiply");
  }
}

/** Map from DiagramIndex.id to DilutedTrace for all random index combinations */
template <size_t rvecs>
using DilutedTraces = std::map<ssize_t, std::vector<DilutedTrace<rvecs>>>;
