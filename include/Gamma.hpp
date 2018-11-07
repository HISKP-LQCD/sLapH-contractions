#pragma once

#include "typedefs.hpp"

#include <vector>

struct gamma_lookup {
  std::array<int, 4> row;
  std::array<Complex, 4> value;
};

std::vector<gamma_lookup> make_gamma();

static std::vector<gamma_lookup> const gamma_vec = make_gamma();
