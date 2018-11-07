#include "dilution-iterator.hpp"

int main() {
  test_dilution_scheme(6, 2, DilutionType::block);
  test_dilution_scheme(6, 2, DilutionType::interlace);

  test_dilution_scheme(48, 4, DilutionType::block);
  test_dilution_scheme(48, 4, DilutionType::interlace);

  test_dilution_scheme(48, 24, DilutionType::block);
}
