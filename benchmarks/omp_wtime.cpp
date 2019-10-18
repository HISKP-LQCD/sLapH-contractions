#include <benchmark/benchmark.h>
#include <omp.h>

static void BM_omp_wtime(benchmark::State &state) {
  volatile double time;

  for (auto _ : state) {
    time __attribute__((unused)) = omp_get_wtime();
  }
}

BENCHMARK(BM_omp_wtime);
