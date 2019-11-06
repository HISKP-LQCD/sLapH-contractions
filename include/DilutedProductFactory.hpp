#pragma once

#include "DilutedFactor.hpp"
#include "DilutedFactorFactory.hpp"

class DilutedProductFactoryQ0Q2 {
 public:
  using TimeKey = std::array<int, 4>;
  using QnKey = std::array<ssize_t, 2>;
  using Value = DilutedFactorsMap<2>;
  using FullKey = std::pair<TimeKey, QnKey>;

  DilutedProductFactoryQ0Q2(DilutedFactorFactory<DilutedFactorType::Q0> &factory_q0,
                            DilutedFactorFactory<DilutedFactorType::Q2> &factory_q2)
      : factory_q0_(factory_q0), factory_q2_(factory_q2) {}

  void request(TimeKey const &time_key, QnKey const &qn_key);

  void build_all();

  std::vector<DilutedFactor> const &get(TimeKey const &time_key, QnKey const &key) {
    return Q0Q2_.at(time_key).at(key);
  }

  void clear() { Q0Q2_.clear(); }

 private:
  void build(TimeKey const &time_key, std::array<ssize_t, 2> const &key);

  std::map<TimeKey, Value> Q0Q2_;

  std::set<FullKey> requests_;

  DilutedFactorFactory<DilutedFactorType::Q0> &factory_q0_;
  DilutedFactorFactory<DilutedFactorType::Q2> &factory_q2_;
};
