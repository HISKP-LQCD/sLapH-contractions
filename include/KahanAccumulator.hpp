#pragma once

#include "ComplexProduct.hpp"

#include <atomic>
#include <cmath>

/**
 * Kahan summation implementation.
 *
 * The algorithm is taken from
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 *
 * @tparam Numeric A numeric type which supports addition and subtraction. The
 * most sensible values are `double` and `std::complex<double>`.
 */
template <typename Numeric>
class KahanAccumulator {
 public:
  Numeric value() const { return sum_ + c_; }

  /**
   * Algorithm taken from
   * https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements.
   */
  KahanAccumulator &operator+=(Numeric right) {
    Numeric const t = sum_ + right;
    if (std::abs(sum_) >= std::abs(right)) {
      c_ += (sum_ - t) + right;
    } else {
      c_ += (right - t) + sum_;
    }
    sum_ = t;
    return *this;
  }

  KahanAccumulator operator+(Numeric const right) const {
    KahanAccumulator left = *this;
    return left += right;
  }

  KahanAccumulator operator+=(KahanAccumulator const right) {
    sum_ += right.sum_;
    c_ += right.c_;
    return *this;
  }

 private:
  Numeric sum_ = 0.0;
  Numeric c_ = 0.0;
};

template <typename Numeric, typename Internal = Numeric>
class NativeAccumulator {
 public:
  NativeAccumulator() : sum_(0.0) {}

  Numeric value() const { return sum_; }

  NativeAccumulator operator+=(Numeric const right) {
    sum_ += right;
    return *this;
  }

  NativeAccumulator operator+(Numeric const right) const {
    NativeAccumulator left = *this;
    return left += right;
  }

  NativeAccumulator operator+=(NativeAccumulator const right) {
    sum_ += right.sum_;
    return *this;
  }

 private:
  Internal sum_;
};

template <typename Numeric, typename Internal = Numeric>
#ifdef SLAPH_KAHAN
using Accumulator = KahanAccumulator<Numeric, Internal>;
#else
using Accumulator = NativeAccumulator<Numeric, Internal>;
#endif

class AtomicAccumulator {
 public:
  AtomicAccumulator() : real_(0.0), imag_(0.0) {}

  Complex value() const { return {real_, imag_}; }

  void operator+=(Complex const right) {
#pragma omp atomic
    real_ += right.real();
#pragma omp atomic
    imag_ += right.imag();
  }

 private:
  double real_;
  double imag_;
};
