#include "DilutedTrace.hpp"

#include "timings.hpp"

#include <iomanip>

Complex inner_product(DilutedTraces const &left_vec, DilutedTraces const &right_vec) {
  TimingScope<4> timing_scope("inner_product(DT, DT)");

  assert(left_vec.traces.size() > 0);
  assert(right_vec.traces.size() > 0);

  int num_summands = 0;
  Accumulator<Complex> result;

  for (auto const &left : left_vec.traces) {
    Accumulator<Complex> right_sum;

    for (auto const &right : right_vec.traces) {
      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if ((left.used_rnd_ids & right.used_rnd_ids) != 0u) {
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

    result += left.data * right_sum.value();
  }

  if (num_summands == 0) {
    throw std::runtime_error(
        "compcomp_t inner_product(vector<DilutedTrace>, "
        "vector<DilutedTrace>) has an empty result. Not enough random vectors?");
  }

  return result.value() / static_cast<double>(num_summands);
}

Complex inner_product(DilutedTraces const &left_vec,
                      DilutedTraces const &middle_vec,
                      DilutedTraces const &right_vec) {
  TimingScope<4> timing_scope("inner_product(DT, DT, DT)");

  assert(left_vec.traces.size() > 0);
  assert(middle_vec.traces.size() > 0);
  assert(right_vec.traces.size() > 0);

  int num_summands = 0;
  Accumulator<Complex> result;

  for (auto const &left : left_vec.traces) {
    Accumulator<Complex> middle_sum;
    for (auto const &middle : middle_vec.traces) {
      Accumulator<Complex> right_sum;

      for (auto const &right : right_vec.traces) {
        if ((left.used_rnd_ids & middle.used_rnd_ids & right.used_rnd_ids) != 0u) {
          continue;
        }

        std::cout
          << left.used_rnd_ids << "\t"
          << std::setprecision(std::numeric_limits<double>::max_digits10)
          << left.data.real() << "\t"
          << std::setprecision(std::numeric_limits<double>::max_digits10)
          << left.data.imag() << "\t"
          << middle.used_rnd_ids << "\t"
          << std::setprecision(std::numeric_limits<double>::max_digits10)
          << middle.data.real() << "\t"
          << std::setprecision(std::numeric_limits<double>::max_digits10)
          << middle.data.imag() << "\t"
          << right.used_rnd_ids << "\t"
          << std::setprecision(std::numeric_limits<double>::max_digits10)
          << right.data.real() << "\t"
          << std::setprecision(std::numeric_limits<double>::max_digits10)
          << right.data.imag() << "\n";

        right_sum += right.data;
        ++num_summands;
      }
      middle_sum += middle.data * right_sum.value();
    }
    result += left.data * middle_sum.value();
  }

  if (num_summands == 0) {
    throw std::runtime_error(
        "compcomp_t inner_product(vector<DilutedTrace>, "
        "vector<DilutedTrace>) has an empty result. Not enough random vectors?");
  }

  return result.value() / static_cast<double>(num_summands);
}

std::vector<DilutedTrace> factor_to_trace(std::vector<DilutedFactor> const &left_vec,
                                          std::vector<DilutedFactor> const &right_vec) {
  TimingScope<4> timing_scope("factor_to_trace(v<DF>, v<DF>)");

  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  std::vector<DilutedTrace> result_vec;

  for (auto const &left : left_vec) {
    TimingScope<5> timing_scope("factor_to_trace(v<DF>, v<DF>): right side");
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

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
      if ((left.used_rnd_ids & right.used_rnd_ids) != 0u) {
        continue;
      }

      /*
      if (((1u << left.ric.first) & right.used_rnd_ids) != 0) {
        continue;
      }
      if (((1u << right.ric.second) & left.used_rnd_ids) != 0) {
        continue;
      }
      */

      // We want to keep track of the indices that have been contracted away. These are
      // all the ones from the left factor, all the ones from the right factor and the one
      // that we are contracting over right now.
      UsedRnd used = left.used_rnd_ids | right.used_rnd_ids | (1u << inner_rnd_id) |
                     (1u << outer_rnd_id);

      // The right sides that we encounter at this point have the same left and right
      // random vector indices. They may differ in the set of used random vector indices.
      // But since we do not plan to contract the result with more DilutedFactor
      // instances, we do not care to preserve the information about the used random
      // vector indices. Therefore we can sum all these elements up to have less
      // multiplications to do.
      result_vec.push_back({DilutedTrace::Data{(left.data * right.data).trace()}, used});
    }
  }  // for(left_vec)

  if (result_vec.size() == 0) {
    throw std::runtime_error(
        "vector<DilutedTrace> factor_to_trace(vector<DilutedFactor>, "
        "vector<DilutedFactor>) has an empty result. Not enough random vectors?");
  }

  return result_vec;
}

std::vector<DilutedTrace> factor_to_trace(std::vector<DilutedFactor> const &vec) {
  TimingScope<4> timing_scope("factor_to_trace(v<DF>)");

  assert(vec.size() > 0);

  std::vector<DilutedTrace> result_vec;

  for (auto const &elem : vec) {
    // We only want to use diagonal elements.
    if (elem.ric.first != elem.ric.second) {
      continue;
    }

    UsedRnd used = elem.used_rnd_ids | (1u << elem.ric.first);
    DilutedTrace result = {elem.data.trace(), used};
    result_vec.push_back(result);
  }

  if (result_vec.size() == 0) {
    throw std::runtime_error(
        "vector<DilutedTrace> factor_to_trace(vector<DilutedFactor>) has an empty "
        "result.");
  }

  return result_vec;
}
