#include <benchmark/benchmark.h>
#include <omp.h>

static void BM_omp_wtime(benchmark::State &state) {
  volatile double time __attribute__((unused));

  for (auto _ : state) {
    time = omp_get_wtime();
  }
}

BENCHMARK(BM_omp_wtime);
