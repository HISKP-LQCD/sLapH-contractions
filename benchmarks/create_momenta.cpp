#include "OperatorsForMesons.hpp"

#include <benchmark/benchmark.h>

static void BM_create_momenta_old(benchmark::State &state) {
  ssize_t const Lx = 32, Ly = 32, Lz = 32;
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup{{0, {1, 2, 3}, {}}};
  array_cd_d2 momentum(boost::extents[1][Lx * Ly * Lz]);

  for (auto _ : state) {
    create_momenta(Lx, Ly, Lz, vdaggerv_lookup, momentum);
  }
}

BENCHMARK(BM_create_momenta_old);

static void BM_create_momenta_new(benchmark::State &state) {
  ssize_t const Lx = 32, Ly = 32, Lz = 32;
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup{{0, {1, 2, 3}, {}}};
  array_cd_d2 momentum(boost::extents[1][Lx * Ly * Lz]);

  for (auto _ : state) {
    create_momenta_new(Lx, Ly, Lz, vdaggerv_lookup, momentum);
  }
}

BENCHMARK(BM_create_momenta_new);
