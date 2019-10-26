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

  void request(TimeKey const &time_key, QnKey const &qn_key) {
    requests_.insert(std::make_pair(time_key, qn_key));
  }

  void build_all() {
    TimingScope<4> timing_scope("DilutedProductFactoryQ0Q2::build_all");

    std::vector<FullKey> unique_requests;
    unique_requests.reserve(requests_.size());
    for (auto const &request : requests_) {
      if (Q0Q2_.count(request.first) == 0 ||
          Q0Q2_.at(request.first).count(request.second) == 0) {
        unique_requests.push_back(request);
        Q0Q2_[request.first];
      }
    }

    for (int i = 0; i < ssize(unique_requests); ++i) {
      auto const &request = unique_requests[i];
      build(request.first, request.second);
    }

    requests_.clear();
  }

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
