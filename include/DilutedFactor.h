/*! @file Operator.h
 *  Declaration of Operator multiplication
 *
 *  @author Markus Werner
 *
 */

#pragma once

#include <iosfwd>
#include <set>
#include <sstream>
#include <vector>

#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "QuarkLineBlock.h"
#include "typedefs.h"
#include "global_data.h"

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

/*! Product yielding the off-diagonal elements.

  From the sets of DilutedFactor elements, the product set of DilutedFactor is build such
  that it only contains elements with _unequal_ left and right random vector index. This
  set is intended to be used as an intermediate result.
  */
template <size_t rvecs1, size_t rvecs2>
std::vector<DilutedFactor<rvecs1 + rvecs2 + 1>> operator*(
    std::vector<DilutedFactor<rvecs1>> const &left_vec,
    std::vector<DilutedFactor<rvecs2>> const &right_vec) {
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

      result_vec.push_back({Eigen::MatrixXcd{left.data * right.data},
                            std::make_pair(left.ric.first, right.ric.second),
                            used});
    }
  }

  return result_vec;
}

template <size_t rvecs>
inline cmplx operator+(DilutedTrace<rvecs> const &df, cmplx const &c) {
  return c + df.data;
}

template <size_t rvecs>
inline cmplx operator+(cmplx const &c, DilutedTrace<rvecs> const &df) {
  return df + c;
}

template <int n, size_t rvecs>
using OperatorToFactorMap =
    std::map<std::array<size_t, n>, std::vector<DilutedFactor<rvecs>>>;

#if 0
template <int n1, int n2>
OperatorToFactorMap<n1 + n2> operator*(OperatorToFactorMap<n1> const &left_map,
                                       OperatorToFactorMap<n2> const &right_map) {
  OperatorToFactorMap<n1 + n2> result;

  for (auto const &left : left_map) {
    for (auto const &right : right_map) {
      // Concatenate the two keys from the left and right element into the new key.
      typename OperatorToFactorMap<n1 + n2>::key_type key;
      auto out_it =
          std::copy(std::begin(left.first), std::end(left.first), std::begin(key));
      std::copy(std::begin(right.first), std::end(right.first), out_it);

      // Do the actual multiplication.
      result[key] = left.second * right.second;
    }
  }
}
#endif

// Proposed:
// template <int n>
// using OperatorToFactorMap = std::map<std::array<QuantumNumbers, n>,
// std::vector<DilutedFactor>>;

template <size_t n, size_t rvecs>
std::string to_string(typename OperatorToFactorMap<n, rvecs>::key_type const &array) {
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
void print(OperatorToFactorMap<n, rvecs> const &otfm) {
  std::cout << "OperatorToFactorMap, size = " << otfm.size() << "\n";
  for (auto const &elem : otfm) {
    std::cout << "  " << to_string<n, rvecs>(elem.first) << " -> "
              << "std::vector(size = " << elem.second.size() << ")\n";
  }
}

/*! @todo Be more restrictive with lookup tables. .Q2V etc. is enough */
template <QuarkLineType qlt>
void check_random_combinations(std::string const &diagram,
                               std::vector<size_t> const &lookup,
                               std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup);

/*! Multiply (Q2V*rVdaggerV) and take trace
 *  - corrC
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
std::vector<cmplx> trace(QuarkLineBlock<qlt1> const &quarkline1,
                         QuarkLineBlock<qlt2> const &quarkline2,
                         int const t1,
                         int const b2,
                         int const t2,
                         std::vector<size_t> const &lookup,
                         std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                         std::vector<size_t> const &ric_ids,
                         int const gamma,
                         size_t const dilE,
                         size_t const dilD);

/*! Multiply (Q1*Q1) and take trace
 *  - corr0
 */
std::vector<cmplx> trace(std::vector<Eigen::MatrixXcd> const &quarkline1,
                         std::vector<Eigen::MatrixXcd> const &quarkline2,
                         std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                         std::vector<size_t> const &ric_ids);

/*! Multiply two traces of two Quarklines each: tr(QQ) * tr(QQ)
 */
compcomp_t trtr(std::vector<cmplx> const &factor1,
                std::vector<cmplx> const &factor2,
                std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                std::vector<size_t> const &ric_ids);

/******************************************************************************/

/*! Create vector<MatrixXcd> with Q1 for all rnd_vecs
 *  - C1
 *  - C3c
 *  - C30
 */
void Q1(std::vector<Eigen::MatrixXcd> &result,
        std::vector<Eigen::MatrixXcd> const &quarklines,
        std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
        std::vector<size_t> const &ric_ids,
        size_t const dilE,
        size_t const dilD);

#if 0
void Q1xQ1(std::vector<DilutedFactor> &result,
           std::vector<DilutedFactor> const &quarkline1,
           std::vector<DilutedFactor> const &quarkline2,
           std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
           std::vector<size_t> const ric_ids,
           size_t const dilE,
           size_t const dilD);
#endif

/*! Create vector<MatrixXcd> with Q0*Q2 for all rnd vecs not equal
 *  - (corrC)
 *  - C4cB
 *  - C4cC
 *  - C3c
 */
void rVdaggerVrxQ2(std::vector<Eigen::MatrixXcd> &result,
                   std::vector<Eigen::MatrixXcd> const &quarkline1,
                   std::vector<Eigen::MatrixXcd> const &quarkline2,
                   std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                   std::vector<size_t> const &ric_ids,
                   size_t const dilE,
                   size_t const dilD);

/*! Multiply (QQ)*(Q) and take trace
 *  - C3c
 *  - C30
 *  @calls M1xM2 for Optimization
 */
cmplx trace_3pt(std::vector<Eigen::MatrixXcd> const &M1,
                std::vector<Eigen::MatrixXcd> const &M2,
                std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                std::vector<size_t> const &ric_ids,
                size_t const dilE,
                size_t const dilD);

/*! Multiply (QQ)*(QQ). or
 *  and take trace
 *  - C4cC
 *  - C4cB
 *  - C40C
 *  - C40B
 *  @calls M1xM2 for Optimization
 */
cmplx trace(std::vector<Eigen::MatrixXcd> const &M1,
            std::vector<Eigen::MatrixXcd> const &M2,
            std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
            std::vector<size_t> const &ric_ids,
            size_t const dilE,
            size_t const dilD);

template <size_t rvecs1, size_t rvecs2>
cmplx trace(std::vector<DilutedFactor<rvecs1>> const &left_vec,
            std::vector<DilutedFactor<rvecs2>> const &right_vec) {
  cmplx result(0.0, 0.0);

  for (auto const &left : left_vec) {
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

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
    }

    auto const &product = left.data * right_sum;
    result += product.trace();
  }

  return result;
}

template <size_t rvecs1, size_t rvecs2>
std::vector<DilutedTrace<rvecs1 + rvecs2 + 2>> factor_to_trace(
    std::vector<DilutedFactor<rvecs1>> const &left_vec,
    std::vector<DilutedFactor<rvecs2>> const &right_vec) {
  int constexpr rvecs_result = rvecs1 + rvecs2 + 2;
  std::vector<DilutedTrace<rvecs_result>> result_vec;

  for (auto const &left : left_vec) {
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

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
  }

  return result_vec;
}

template <size_t rvecs>
std::vector<DilutedTrace<rvecs + 1>> factor_to_trace(
    std::vector<DilutedFactor<rvecs>> const &vec) {
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

  return result_vec;

}

template <size_t rvecs1, size_t rvecs2>
compcomp_t inner_product(std::vector<DilutedTrace<rvecs1>> const &left_vec,
                         std::vector<DilutedTrace<rvecs2>> const &right_vec) {
  //! @TODO Pull out this magic number.
  auto constexpr rnd_vec_count = 5;

  compcomp_t result(0.0, 0.0, 0.0, 0.0);

  for (auto const &left : left_vec) {
    cmplx right_sum(0.0, 0.0);

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
    }

    result.rere += left.data.real() * right_sum.real();
    result.reim += left.data.real() * right_sum.imag();
    result.imre += left.data.imag() * right_sum.real();
    result.imim += left.data.imag() * right_sum.imag();
  }

  return result;
}

int constexpr max_flavor = 8;
using UpTo = boost::container::static_vector<RndId, max_flavor>;

template <size_t rvecs>
UpTo get_max_used(DilutedFactor<rvecs> const &df, UpTo const &rnd_offset) {
  SmallVectorRndId<rvecs + 2> used = df.used_rnd_ids;
  merge_push_back(used, df.ric.first);
  merge_push_back(used, df.ric.second);

  UpTo result(rnd_offset.size() - 1, 0);

  // Iterate through every random vector index that is either used internally or
  // externally.
  for (auto const id : used) {
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

template <typename T, size_t rvecs1, size_t rvecs2>
T by_rnd_vec_count(std::function<T(std::vector<DilutedFactor<rvecs1>> const &,
                                   std::vector<DilutedFactor<rvecs2>> const &)> func,
                   std::vector<DilutedFactor<rvecs1>> const &left_vec,
                   std::vector<DilutedFactor<rvecs2>> const &right_vec) {
  // Extract the number of random vectors per quark flavor.
  auto const &quarks = GlobalData::Instance()->quarks();
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

  std::map<UpTo, T> results;

  // We start with the case where in every flavor we only use 1 random vector.
  UpTo upto(quarks.size(), 1);

  // We want to go until we have used all available random vectors.
  UpTo upto_max = rnd_count;

  while (upto != upto_max) {
    // We only want to take those diluted factors that make use of the highest allowed
    // random vector index for each flavor, but nothing more and nothing less. This avoids
    // repeated use of the same `DilutedFactor`.
    std::vector<DilutedFactor<rvecs1>> left_filtered;
    std::copy_if(std::begin(left_vec),
                 std::end(left_vec),
                 std::back_inserter(left_filtered),
                 [&rnd_offset, &upto](DilutedFactor<rvecs1> const &df) {
                   return get_max_used(df, rnd_offset) == upto;
                 });

    // Same for the right.
    std::vector<DilutedFactor<rvecs2>> right_filtered;
    std::copy_if(std::begin(right_vec),
                 std::end(right_vec),
                 std::back_inserter(right_filtered),
                 [&rnd_offset, &upto](DilutedFactor<rvecs2> const &df) {
                   return get_max_used(df, rnd_offset) == upto;
                 });

    // Now we apply the function to the limited set.
    results[upto] = func(left_filtered, right_filtered);

    // Iterate to the next element. Do this by increasing the lowest digit. If that
    // overflows, reset it and increase the next digit by one.
    ++upto[0];
    for (int i = 1; i < upto.size(); ++i) {
      if (upto[i - 1] == rnd_count[i - 1] + 1) {
        upto[i - 1] = 1;
        ++upto[i];
      }
    }
  }

  return results;
}
