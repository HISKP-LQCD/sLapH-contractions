#include "dilution-iterator.hpp"

#include <iostream>
#include <vector>

int get_time_delta(BlockIterator const &slice_pair, int const Lt) {
  return abs((slice_pair.sink() - slice_pair.source() - Lt) % Lt);
}

int main() {
  auto const Lt = 96;
  auto const dilT = 3;
  auto const time_slice_divisor = 3;
  auto const time_slice_remainder = 0;

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  std::vector<int> counts(Lt, 0);

  for (int b = 0; b < dilution_scheme.size(); ++b) {
    auto const block_pair = dilution_scheme[b];
    for (auto const slice_pair : block_pair) {
      int const t = get_time_delta(slice_pair, Lt);
      if (slice_pair.source() % time_slice_divisor != time_slice_remainder) {
        continue;
      }

      ++counts[t];
    }
  }

  for (auto t = 0u; t < counts.size(); ++t) {
    std::cout << t << "\t" << counts[t] << "\n";
  }
}
