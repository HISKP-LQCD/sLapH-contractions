#include "DilutedProductFactory.hpp"

#include "timings.hpp"

void DilutedProductFactoryQ0Q2::request(TimeKey const &time_key, QnKey const &qn_key) {
  TimingScope<5> timing_scope("DilutedProductFactoryQ0Q2::request");

  // Extract time keys.
  int constexpr nt1 = 3;
  int constexpr nt2 = 1;
  std::array<int, nt1> time_key1;
  std::array<int, nt2> time_key2;
  std::copy_n(std::begin(time_key) + 0, nt1, std::begin(time_key1));
  std::copy_n(std::begin(time_key) + nt1, nt2, std::begin(time_key2));

  factory_q2_.request(time_key1);
  factory_q0_.request(time_key2);

  requests_.insert(std::make_pair(time_key, qn_key));
}

void DilutedProductFactoryQ0Q2::build_all() {
  TimingScope<3> timing_scope("DilutedProductFactoryQ0Q2::build_all");

  std::vector<FullKey> unique_requests;
  unique_requests.reserve(requests_.size());
  for (auto const &request : requests_) {
    if (Q0Q2_.count(request.first) == 0 ||
        Q0Q2_.at(request.first).count(request.second) == 0) {
      unique_requests.push_back(request);
      Q0Q2_[request.first][request.second];
    }
  }

  requests_.clear();

  //#pragma omp parallel for
  for (int i = 0; i < ssize(unique_requests); ++i) {
    auto const &request = unique_requests[i];
    build(request.first, request.second);
    assert(Q0Q2_[request.first][request.second].size() > 0);
  }
}

template <typename Type, size_t size>
std::ostream &operator<<(std::ostream &os, std::array<Type, size> const &array) {
  os << "[";
  bool first = true;
  for (auto const &elem : array) {
    if (!first) {
      os << ", ";
    }
    os << elem;
    first = false;
  }
  os << "]";

  return os;
}

void DilutedProductFactoryQ0Q2::build(TimeKey const &time_key, QnKey const &key) {
  TimingScope<5> timing_scope("DilutedProductFactoryQ0Q2::build");

  // Extract time keys.
  int constexpr nt1 = 3;
  int constexpr nt2 = 1;
  std::array<int, nt1> time_key1;
  std::array<int, nt2> time_key2;
  std::copy_n(std::begin(time_key) + 0, nt1, std::begin(time_key1));
  std::copy_n(std::begin(time_key) + nt1, nt2, std::begin(time_key2));

  multiply<1, 1>(Q0Q2_[time_key], key, factory_q2_[time_key1], factory_q0_[time_key2]);

  if (Q0Q2_[time_key][key].size() == 0) {
#pragma omp critical(cerr)
    std::cerr << "Failed Q0Q2: " << time_key << " " << key << std::endl;
    abort();
  }
}
