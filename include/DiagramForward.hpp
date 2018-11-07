#pragma once

#include "typedefs.hpp"

template <typename Numeric>
class DiagramNumeric;

template <typename Numeric>
struct DiagramTraits {
  using Diagram = DiagramNumeric<Numeric>;
};

class DiagramParts;
