#pragma once

/**
 * Cartesian product similar to Python's `itertools.product`.
 *
 * Give it a list of lengths and it will give you a sequence with a set of
 * indices. This is different from the Python implementation which directly
 * works with the data. We use indices here because of all that fancy lookup
 * stuff.
 *
 * Damn, there is not a single occurence of “lookup” in this class. So let me
 * try again: Give this class the lenths of your lookup tables and it will
 * return a sequence of lookup arrays which give you the indices which you can
 * use to look up the correct permutation of your data in your lookup data
 * structures.
 *
 * There is the Boost HANA library which also provides a cartesian product. It
 * is much nicer than this one here because it works with the actual data and
 * is resolved during compile time. The problem is that we want it during
 * runtime here, so we cannot use this. Also it requires C++14, which we cannot
 * have here as this is a _production code_ and not some playground. You never
 * know what ancient compiler needs to be supported, so let us be grateful for
 * having C++11.
 */
class CartesianProduct {
 public:
  using Lengths = std::vector<ssize_t>;

  class Iterator {
   public:
    Iterator(Lengths const &lengths, int i) : lengths_(lengths), i_(i) {}

    std::vector<size_t> operator*() {
      std::vector<size_t> combination(lengths_.size());
      auto ii = i_;
      for (size_t j = 0; j != lengths_.size(); ++j) {
        combination[j] = ii % lengths_.at(j);
        ii /= lengths_.at(j);
      }
      return combination;
    }

    Iterator operator++() {
      ++i_;
      return *this;
    }

    bool operator!=(Iterator const &other) {
      return lengths_ != other.lengths_ || i_ != other.i_;
    }

   private:
    Lengths const &lengths_;
    ssize_t i_;
  };

  CartesianProduct(Lengths const &lengths) : lengths_(lengths) {}

  ssize_t size() {
    ssize_t combinations = 1;
    for (auto const length : lengths_) {
      combinations *= length;
    }
    return combinations;
  }

  Iterator begin() { return Iterator(lengths_, 0); }

  Iterator end() { return Iterator(lengths_, size()); }

 private:
  Lengths const lengths_;
};
