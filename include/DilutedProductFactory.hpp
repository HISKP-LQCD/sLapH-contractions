#pragma once

#include "DilutedFactor.hpp"
#include "DilutedFactorFactory.hpp"

class DilutedProductFactoryQ0Q2 {
 public:
  using Key = std::array<int, 4>;
  using Value = DilutedFactorsMap<2>;

  DilutedProductFactoryQ0Q2(DilutedFactorFactory<DilutedFactorType::Q0> &factory_q0,
                            DilutedFactorFactory<DilutedFactorType::Q2> &factory_q2)
      : factory_q0_(factory_q0), factory_q2_(factory_q2) {}

  std::vector<DilutedFactor> const &get(Key const &time_key,
                                        std::array<ssize_t, 2> const &key) {
    TimingScope<4> timing_scope("DilutedProductFactoryQ0Q2::get");

      if (Q0Q2_.count(time_key) == 0 || Q0Q2_.at(time_key).count(key) == 0) {
        build(time_key, key);
      }

    return Q0Q2_.at(time_key).at(key);
  }

  void clear() { Q0Q2_.clear(); }

  void build(Key const &time_key, std::array<ssize_t, 2> const &key);

 private:
  std::map<Key, Value> Q0Q2_;

  DilutedFactorFactory<DilutedFactorType::Q0> &factory_q0_;
  DilutedFactorFactory<DilutedFactorType::Q2> &factory_q2_;
};
