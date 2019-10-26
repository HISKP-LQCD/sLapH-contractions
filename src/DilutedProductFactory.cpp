#include "DilutedProductFactory.hpp"

#include "timings.hpp"

void DilutedProductFactoryQ0Q2::build(Key const &time_key,
                                      std::array<ssize_t, 2> const &key) {
  TimingScope<5> timing_scope("DilutedProductFactoryQ0Q2::build");

  // Extract time keys.
  int constexpr nt1 = 3;
  int constexpr nt2 = 1;
  std::array<int, nt1> time_key1;
  std::array<int, nt2> time_key2;
  std::copy_n(std::begin(time_key) + 0, nt1, std::begin(time_key1));
  std::copy_n(std::begin(time_key) + nt1, nt2, std::begin(time_key2));

  factory_q2_.request(time_key1);
  factory_q0_.request(time_key2);
  factory_q2_.build_all();
  factory_q0_.build_all();

  multiply<1, 1>(Q0Q2_[time_key], key, factory_q2_[time_key1], factory_q0_[time_key2]);
}
