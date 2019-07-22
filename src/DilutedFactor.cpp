#include "DilutedFactor.hpp"

#include "KahanAccumulator.hpp"
#include "timings.hpp"

std::vector<DilutedFactor> operator*(std::vector<DilutedFactor> const &left_vec,
                                     std::vector<DilutedFactor> const &right_vec) {
  TimingScope<4> timing_scope("operator*(v<DF>, v<DF>");

  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  std::vector<DilutedFactor> result_vec;

  for (auto const &left : left_vec) {
    auto const inner_rnd_id = left.ric.second;

    for (auto const &right : right_vec) {
      // The inner indices need to match because the product would not make
      // sense otherwise.
      bool const is_allowed =
          inner_rnd_id == right.ric.first && left.ric.first != right.ric.second;
      if (!is_allowed) {
        continue;
      }

      if (left.ric.first == left.ric.second) {
        continue;
      }

      if (right.ric.first == right.ric.second) {
        continue;
      }

      if (((1u << left.ric.first) & right.used_rnd_ids) != 0) {
        continue;
      }
      if (((1u << right.ric.second) & left.used_rnd_ids) != 0) {
        continue;
      }

      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if ((left.used_rnd_ids & right.used_rnd_ids) != 0u) {
        continue;
      }

      // We want to keep track of the indices that have been contracted away. These are
      // all the ones from the left factor, all the ones from the right factor and the one
      // that we are contracting over right now.
      UsedRnd used = left.used_rnd_ids | right.used_rnd_ids | (1u << inner_rnd_id);

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

Complex trace(std::vector<DilutedFactor> const &left_vec,
              std::vector<DilutedFactor> const &right_vec) {
  TimingScope<4> timing_scope("trace(v<DF>, v<DF>)");

  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  int num_summands = 0;
  Accumulator<Complex> result;

  for (auto const &left : left_vec) {
    TimingScope<4> timing_scope("trace(v<DF>, v<DF>): right sum");

    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

    // XXX (MU 2019-01-03): This variable is an accumulator, but it runs with
    // double precision in all cases. The discussion has been started at
    // https://github.com/HISKP-LQCD/sLapH-contractions/pull/87#issuecomment-451157721.
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
      if ((left.used_rnd_ids & right.used_rnd_ids) == 0u) {
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

    auto const &product = left.data * right_sum;
    result += product.trace();

  }  // for(left_vec)

  if (num_summands == 0) {
    throw std::runtime_error(
        "cmplx trace(vector<DilutedFactor>, vector<DilutedFactor>) has used zero "
        "summands.");
  }

  return result.value() / static_cast<double>(num_summands);
}

std::ostream &operator<<(std::ostream &os, DilutedFactor const &df) {
  os << "DilutedFactor{{";
  auto const &m = df.data;
  for (int row = 0; row < m.rows(); ++row) {
    if (row != 0) {
      os << "; ";
    }
    for (int col = 0; col < m.cols(); ++col) {
      if (col != 0) {
        os << ", ";
      }
      os << m(row, col);
    }
  }
  os << "}, {" << static_cast<int>(df.ric.first) << ", " << static_cast<int>(df.ric.second)
     << "}, {";
  auto used = df.used_rnd_ids;
  int base = 0;
  bool output = false;
  while (used != 0) {
    if (used % 2 == 1) {
    if (output) {
      os << ", ";
    }
      os << base;
      output = true;
    }
    used /= 2;
    ++base;
  }
  os << "}}";
  return os;
}
