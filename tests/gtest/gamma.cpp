#include "Gamma.hpp"

#include <gtest/gtest.h>

namespace {
std::complex<double> const I{0.0, 1.0};
}

Eigen::Matrix4cd sparse_to_eigen_dense(gamma_lookup const sparse) {
  Eigen::Matrix4cd dense;
  for (ssize_t col = 0; col < 4; ++col) {
    dense(sparse.row[col], col) = sparse.value[col];
  }
  return dense;
}

// TODO: How to best reuse tmlqcd_gamma matrices? struct?
TEST(gamma, checkBasis) {
  // Gamma matrices as defined in tmLQCD
  std::array<Eigen::Matrix4cd, 4> tmlqcd_gamma;
  tmlqcd_gamma[0] <<  0,  0,  1,  0,
                      0,  0,  0,  1,
                      1,  0,  0,  0,
                      0,  1,  0,  0;

  tmlqcd_gamma[1] <<  0,  0,  0,  I,
                      0,  0,  I,  0,
                      0, -I,  0,  0,
                     -I,  0,  0,  0;

  tmlqcd_gamma[2] <<  0,  0,  0,  1,
                      0,  0, -1,  0,
                      0, -1,  0,  0,
                      1,  0,  0,  0;

  tmlqcd_gamma[3] <<  0,  0,  I,  0,
                      0,  0,  0, -I,
                     -I,  0,  0,  0,
                      0,  I,  0,  0;
  //Check basis matrices
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[0]), tmlqcd_gamma[0]);
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[1]), tmlqcd_gamma[1]);
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[2]), tmlqcd_gamma[2]);
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[3]), tmlqcd_gamma[3]);
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[4]), Eigen::Matrix4cd::Identity());
}

TEST(gamma, checkProducts) {
  // Gamma matrices as defined in tmLQCD
  std::array<Eigen::Matrix4cd, 4> tmlqcd_gamma;
  tmlqcd_gamma[0] <<  0,  0,  1,  0,
                      0,  0,  0,  1,
                      1,  0,  0,  0,
                      0,  1,  0,  0;

  tmlqcd_gamma[1] <<  0,  0,  0,  I,
                      0,  0,  I,  0,
                      0, -I,  0,  0,
                     -I,  0,  0,  0;

  tmlqcd_gamma[2] <<  0,  0,  0,  1,
                      0,  0, -1,  0,
                      0, -1,  0,  0,
                      1,  0,  0,  0;

  tmlqcd_gamma[3] <<  0,  0,  I,  0,
                      0,  0,  0, -I,
                     -I,  0,  0,  0,
                      0,  I,  0,  0;

  // Check gamma5
  Eigen::Matrix4cd test_gamma5;
  test_gamma5 = Eigen::Matrix4cd::Identity();
  for (auto const &g : tmlqcd_gamma)
    test_gamma5 *= g;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[5]), test_gamma5);

  // Check linear combinations
  // gamma_0 * gamma_5
  Eigen::Matrix4cd gamma05 = tmlqcd_gamma[0] * test_gamma5;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[6]), gamma05);

  // gamma_1 * gamma_5
  Eigen::Matrix4cd gamma15 = tmlqcd_gamma[1] * test_gamma5;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[7]), gamma15);

  // gamma_2 * gamma_5
  Eigen::Matrix4cd gamma25 = tmlqcd_gamma[2] * test_gamma5;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[8]), gamma25);

  // gamma_3 * gamma_5
  Eigen::Matrix4cd gamma35 = tmlqcd_gamma[3] * test_gamma5;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[9]), gamma35);

  // gamma_0 * gamma_1
  Eigen::Matrix4cd gamma01 = tmlqcd_gamma[0] * tmlqcd_gamma[1];
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[10]), gamma01);

  // gamma_0 * gamma_2
  Eigen::Matrix4cd gamma02 = tmlqcd_gamma[0] * tmlqcd_gamma[2];
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[11]), gamma02);

  // gamma_0 * gamma_3
  Eigen::Matrix4cd gamma03 = tmlqcd_gamma[0] * tmlqcd_gamma[3];
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[12]), gamma03);

  // gamma_0 * gamma_1 * gamma_5
  Eigen::Matrix4cd gamma015 = tmlqcd_gamma[0] * tmlqcd_gamma[1] * test_gamma5;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[13]), gamma015);

  // gamma_0 * gamma_2 * gamma_5
  Eigen::Matrix4cd gamma025 = tmlqcd_gamma[0] * tmlqcd_gamma[2] * test_gamma5;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[14]), gamma025);

  // gamma_0 * gamma_3 * gamma_5
  Eigen::Matrix4cd gamma035 = tmlqcd_gamma[0] * tmlqcd_gamma[3] * test_gamma5;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[15]), gamma035);
}
