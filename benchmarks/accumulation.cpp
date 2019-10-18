#include "KahanAccumulator.hpp"

#include <benchmark/benchmark.h>

#include <random>
#include <vector>

static std::vector<double> make_numbers() {
  std::vector<double> numbers(500);

  std::default_random_engine engine(0);
  std::uniform_real_distribution<double> dist(0, 1);

  for (auto &number : numbers) {
    number = dist(engine);
  }

  return numbers;
}

template<typename TestAccumulator>
static void bm_accumulate(benchmark::State &state) {
  auto const &numbers = make_numbers();

  for (auto _ : state) {
    TestAccumulator acc;
    for (auto const number : numbers) {
      acc += number;
    }
    volatile double value __attribute__((unused)) = acc.value();
  }
}

static void BM_kahan_double(benchmark::State &state) {
  bm_accumulate<KahanAccumulator<double>>(state);
}

static void BM_native_double(benchmark::State &state) {
  bm_accumulate<NativeAccumulator<double>>(state);
}

static void BM_native_float128(benchmark::State &state) {
  bm_accumulate<NativeAccumulator<__float128>>(state);
}

static void BM_native_quad(benchmark::State &state) {
  bm_accumulate<NativeAccumulator<long double>>(state);
}

BENCHMARK(BM_native_quad);
BENCHMARK(BM_native_double);
BENCHMARK(BM_native_float128);
BENCHMARK(BM_kahan_double);
