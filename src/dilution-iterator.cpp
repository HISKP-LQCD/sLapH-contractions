#include "dilution-iterator.hpp"

std::ostream &operator<<(std::ostream &os, DilutionIterator const &di) {
  os << "DilutionIterator(source=" << di.source() << ", sink=" << di.sink() << ")";
  return os;
}

std::ostream &operator<<(std::ostream &os, BlockIterator const &bi) {
  os << "BlockIterator(source=" << bi.source() << ", sink=" << bi.sink()
     << ", source_block=" << bi.source_block() << ", sink_block=" << bi.sink_block()
     << ")";
  return os;
}

void test_dilution_scheme(int const num_slice,
                          int const num_block,
                          DilutionType const type) {
  auto const name = dilution_names.at(type);
  auto const block_size = num_slice / num_block;
  std::cout << "T = " << num_slice << ", T" << name << num_block << " (Morningstar), T"
            << name << block_size << " (Other):\n\n";

  DilutionScheme dilution_scheme(num_slice, block_size, type);
  for (int b = 0; b < dilution_scheme.size(); ++b) {
    auto const blocks = dilution_scheme[b];
    std::cout << std::setw(2) << blocks.source() << " => " << std::setw(2)
              << blocks.sink() << "\n";

    for (auto const slices : blocks) {
      std::cout << "  " << std::setw(2) << slices.source() << " -> " << std::setw(2)
                << slices.sink() << "\n";
    }
  }

  std::cout << "\n\n";
}
