#include "OperatorsForMesons.hpp"

#include <gtest/gtest.h>

TEST(createMomenta, aniso) {
  ssize_t const Lx = 14, Ly = 26, Lz = 34;

  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup{{0, {2, 3, 5}, {}}};

  array_cd_d2 momentum_new(boost::extents[1][Lx * Ly * Lz]);
  array_cd_d2 momentum_old(boost::extents[1][Lx * Ly * Lz]);

  create_momenta(Lx, Ly, Lz, vdaggerv_lookup, momentum_old);
  create_momenta_new(Lx, Ly, Lz, vdaggerv_lookup, momentum_new);

  for (ssize_t i = 0; i < Lx * Ly * Lz; ++i) {
    EXPECT_NEAR(momentum_old[0][i].real(), momentum_new[0][i].real(), 1e-14);
    EXPECT_NEAR(momentum_old[0][i].imag(), momentum_new[0][i].imag(), 1e-14);
  }
}
