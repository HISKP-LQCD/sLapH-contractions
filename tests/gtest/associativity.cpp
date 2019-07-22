#include "DilutedFactor.hpp"

#include <gtest/gtest.h>

#include <iomanip>
#include <random>

static std::vector<DilutedFactor> make_random_diluted_factor(int seed, ssize_t size = 50) {
  std::vector<DilutedFactor> df;
  df.reserve(size);

  std::default_random_engine engine(seed);
  std::uniform_real_distribution<double> real_dist(0, 1);
  std::uniform_int_distribution<int> int_dist(0, 4);

  for (int i = 0; i < size; ++i) {
    Eigen::MatrixXcd m(2, 2);
    m << real_dist(engine), real_dist(engine), real_dist(engine), real_dist(engine);
    DilutedFactor f{m, {int_dist(engine), int_dist(engine)}, {}};
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
  }
    return true;
}

static bool operator<(DilutedFactor const &df1, DilutedFactor const &df2) {
  return df1.ric.first < df2.ric.first || df1.ric.second < df2.ric.second ||
         df1.used_rnd_ids < df2.used_rnd_ids || df1.data < df2.data;
}

static bool operator==(DilutedFactor const &df1, DilutedFactor const &df2) {
  return df1.data == df2.data && df1.ric.first == df2.ric.first &&
         df1.ric.second == df2.ric.second && df1.used_rnd_ids == df2.used_rnd_ids;
}

TEST(DilutedFactor, associativity) {
  auto df1 = make_random_diluted_factor(1);
  auto df2 = make_random_diluted_factor(2);
  auto df3 = make_random_diluted_factor(3);
  auto df4 = make_random_diluted_factor(4);

  //std::sort(std::begin(df1), std::end(df1));
  //std::sort(std::begin(df2), std::end(df2));
  //std::sort(std::begin(df3), std::end(df3));
  //std::sort(std::begin(df4), std::end(df4));

  EXPECT_NE(df1, df2);
  EXPECT_NE(df1, df3);
  EXPECT_NE(df2, df3);

  EXPECT_EQ(ssize(df1), ssize(df2));
  EXPECT_EQ(ssize(df1), ssize(df3));

  auto prod_12_3 = (df1 * df2) * df3;
  auto prod_1_23 = df1 * (df2 * df3);

  //std::sort(std::begin(prod_12_3), std::end(prod_12_3));
  //std::sort(std::begin(prod_1_23), std::end(prod_1_23));

  //if (ssize(prod_12_3) != ssize(prod_1_23)) {
    std::cout << "df1:\n";
    for (auto const &elem : df1) {
      std::cout << elem << "\n";
    }
    std::cout << "df2:\n";
    for (auto const &elem : df2) {
      std::cout << elem << "\n";
    }
    std::cout << "df3:\n";
    for (auto const &elem : df3) {
      std::cout << elem << "\n";
    }
    std::cout << "1 * 2:\n";
    for (auto const &elem : df1 * df2) {
      std::cout << elem << "\n";
    }
    std::cout << "2 * 3:\n";
    for (auto const &elem : df2 * df3) {
      std::cout << elem << "\n";
    }
    std::cout << "(1 * 2) * 3\n";
    for (auto const &elem : prod_12_3) {
      std::cout << elem << "\n";
    }
    std::cout << "1 * (2 * 3)\n";
    for (auto const &elem : prod_1_23) {
      std::cout << elem << "\n";
    }
  //}

  ASSERT_EQ(ssize(prod_12_3), ssize(prod_1_23));

  for (int i = 0; i < ssize(prod_1_23); ++i) {
    auto const &a = prod_1_23[i];
    auto const &b = prod_12_3[i];

    auto const &m1 = a.data;
    auto const &m2 = b.data;

    for (int row = 0; row < m1.rows(); ++row) {
      for (int col = 0; col < m1.cols(); ++col) {
        auto const &c1 = m1(row, col);
        auto const &c2 = m2(row, col);
        EXPECT_NEAR(c1.real(), c2.real(), 1e-10);
        EXPECT_NEAR(c1.imag(), c2.imag(), 1e-10);
      }
    }

    EXPECT_EQ(a.ric.first, b.ric.first);
    EXPECT_EQ(a.ric.second, b.ric.second);
    EXPECT_EQ(a.used_rnd_ids, b.used_rnd_ids);
  }
}
