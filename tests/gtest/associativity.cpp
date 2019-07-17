#include "DilutedFactor.hpp"

#include <gtest/gtest.h>

#include <iomanip>
#include <random>

static std::vector<DilutedFactor> make_random_diluted_factor(int seed, ssize_t size = 5) {
  std::vector<DilutedFactor> df;
  df.reserve(size);

  std::default_random_engine engine(seed);
  std::uniform_real_distribution<double> real_dist(0, 1);
  std::uniform_int_distribution<int> int_dist(0, 4);

  for (int i = 0; i < size; ++i) {
    Eigen::MatrixXcd m(2, 2);
    m << real_dist(engine), real_dist(engine), real_dist(engine), real_dist(engine);
    DilutedFactor f{m, {1 << int_dist(engine), 1 << int_dist(engine)}, {}};
    df.push_back(f);
  }

  return df;
}

static bool operator<(Eigen::MatrixXcd const &m1, Eigen::MatrixXcd const &m2) {
  assert(m1.rows() == m2.rows());
  assert(m2.cols() == m2.cols());

  for (int row = 0; row < m1.rows(); ++row) {
    for (int col = 0; col < m1.cols(); ++col) {
      if (m1(row, col).real() >= m2(row, col).real()) {
        return false;
      }
      if (m1(row, col).imag() >= m2(row, col).imag()) {
        return false;
      }
    }
    return true;
  }
}

static bool operator<(DilutedFactor const &df1, DilutedFactor const &df2) {
  return df1.data < df2.data && df1.ric < df2.ric && df1.used_rnd_ids < df2.used_rnd_ids;
}

static bool operator==(DilutedFactor const &df1, DilutedFactor const &df2) {
  return df1.data == df2.data && df1.ric == df2.ric && df1.used_rnd_ids == df2.used_rnd_ids;
}

static std::ostream &operator<<(std::ostream &os, DilutedFactor const &df) {
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
    if (output) {
      os << ", ";
    }
    if (used % 2 == 1) {
      os << (1 << base);
      output = true;
    }
    used /= 2;
    ++base;
  }
  os << "}}";
  return os;
}

TEST(DilutedFactor, associativity) {
  auto df1 = make_random_diluted_factor(1);
  auto df2 = make_random_diluted_factor(2);
  auto df3 = make_random_diluted_factor(3);

  EXPECT_NE(df1, df2);
  EXPECT_NE(df1, df3);
  EXPECT_NE(df2, df3);

  EXPECT_EQ(ssize(df1), ssize(df2));
  EXPECT_EQ(ssize(df1), ssize(df3));

  auto prod_12_3 = (df1 * df2) * df3;
  auto prod_1_23 = df1 * (df2 * df3);

  EXPECT_EQ(ssize(prod_12_3), ssize(prod_1_23));

  std::sort(std::begin(prod_12_3), std::end(prod_12_3));
  std::sort(std::begin(prod_1_23), std::end(prod_1_23));

  EXPECT_EQ(prod_12_3, prod_1_23);
}
