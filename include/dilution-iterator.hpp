#pragma once

#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <type_traits>

enum class DilutionType { block, interlace };

std::map<DilutionType const, std::string const> const dilution_names = {
    {DilutionType::block, "B"}, {DilutionType::interlace, "I"}};

class DilutionIterator;
class BlockIterator;

std::ostream &operator<<(std::ostream &os, DilutionIterator const &di);
std::ostream &operator<<(std::ostream &os, BlockIterator const &bi);

class BlockIterator {
 public:
  BlockIterator(int const slice_source,
                int const slice_sink,
                int const block_source,
                int const block_sink,
                int const num_slice,
                int const num_block,
                int const pass,
                DilutionType const type,
                bool const one_sink_slice = false)
      : slice_source_(slice_source),
        slice_sink_(slice_sink),
        block_source_(block_source),
        block_sink_(block_sink),
        num_slice_(num_slice),
        num_block_(num_block),
        pass_(pass),
        type_(type),
        one_sink_slice_(one_sink_slice) {}

  BlockIterator operator*() const {
#ifdef DILUTION_ITERATOR_PRINT
#pragma omp critical(cout)
    std::cout << *this << std::endl;
#endif
    return *this;
  }

  BlockIterator operator++() {
    if (type_ == DilutionType::block) {
      int const block_size = num_slice_ / num_block_;
      int const block_source_begin = block_source_ * block_size;
      int const block_source_end = (block_source_ + 1) * block_size;
      int const block_sink_begin = block_sink_ * block_size;
      int const block_sink_end = (block_sink_ + 1) * block_size;

      // Go to the next sink slice.
      if (one_sink_slice_) {
        slice_sink_ += block_size;
      } else {
        ++slice_sink_;
      }

      // If we iterated through all sinks in this block, we need to go to the
      // next source.
      if (slice_sink_ == block_sink_end) {
        slice_sink_ = block_sink_begin;
        ++slice_source_;
      }

      // If we iterated through all the sources in our block, we need to
      // exchange
      // source and sink block and repeat.
      if (slice_source_ == block_source_end) {
        slice_source_ = block_sink_begin;
        slice_sink_ = block_source_begin;
        ++pass_;

        if (block_source_ == block_sink_) {
          ++pass_;
        }

        std::swap(block_source_, block_sink_);
      }
    } else if (type_ == DilutionType::interlace) {
      int const block_size = num_slice_ / num_block_;
      int const block_source_begin = block_source_;
      int const block_sink_begin = block_sink_;

      // Go to the next sink slice.
      if (one_sink_slice_) {
        slice_sink_ += num_slice_;
      } else {
        slice_sink_ += block_size;
      }

      // If we iterated through all sinks in this block, we need to go to the
      // next source.
      if (slice_sink_ >= num_slice_) {
        slice_sink_ = block_sink_begin;
        slice_source_ += block_size;
      }

      // If we iterated through all the sources in our block, we need to
      // exchange
      // source and sink block and repeat.
      if (slice_source_ >= num_slice_) {
        slice_source_ = block_sink_begin;
        slice_sink_ = block_source_begin;
        ++pass_;

        if (block_source_ == block_sink_) {
          ++pass_;
        }

        std::swap(block_source_, block_sink_);
      }
    } else {
      throw std::domain_error("This dilution scheme is not implemented.");
    }

    return *this;
  }

  bool operator!=(BlockIterator const &other) const {
    return block_source_ != other.block_source_ || block_sink_ != other.block_sink_ ||
           num_slice_ != other.num_slice_ || num_block_ != other.num_block_ ||
           pass_ != other.pass_ || type_ != other.type_;
  }

  int source() const { return slice_source_; }

  int sink() const { return slice_sink_; }

  int source_block() const { return block_source_; }

  int sink_block() const { return block_sink_; }

  /**
    The following quark line index is a contraption that assigns each time
    slice within the current block a number in the interval `[0, 2 * dilT)`.
    The entries `[0, dilT)` correspond to the “block → sink“ step and the
    entries `[dilT, 2 * dilT)` corrspond to the other direction, say “sink ←
    block”.

    The following is the code that had been used previously to compute this
    index, only for the block-dilution scheme:

    ```{.cpp}
    if (t1_i == t2_i) {
      if (t1_min != 0)
        id_Q2L_1 = t1 % t1_min;
      else {
        id_Q2L_1 = t1;
      }
    } else {
      if (t1_min != 0)
        id_Q2L_1 = (dir)*dilT + t1 % t1_min;
      else
        id_Q2L_1 = ((dir)*dilT + t1);
    }
    ```

    Now we can use the container/iterator structure for the blocks and
    slices. In case we are in the first part, the source block has not been
    switched with the sink block.
    */
  int qline_id() const {
    auto const block_size = num_slice_ / num_block_;

    if (type_ == DilutionType::block) {
      return slice_source_ % block_size + pass_ * block_size;
    } else if (type_ == DilutionType::interlace) {
      return slice_source_ / block_size + pass_ * block_size;
    } else {
      throw std::domain_error("This dilution scheme is not implemented.");
    }
  }

 private:
  int slice_source_;
  int slice_sink_;
  int block_source_;
  int block_sink_;
  int num_slice_;
  int num_block_;
  int pass_;
  DilutionType type_;
  bool one_sink_slice_;
};

class DilutionIterator {
 public:
  struct Element {
    int source;
    int sink;
  };

  DilutionIterator(int const block_source,
                   int const block_sink,
                   int const num_slice,
                   int const num_block,
                   DilutionType const type)
      : block_source_(block_source),
        block_sink_(block_sink),
        num_slice_(num_slice),
        num_block_(num_block),
        type_(type) {}

  DilutionIterator operator*() const {
#ifdef DILUTION_ITERATOR_PRINT
#pragma omp critical(cout)
    std::cout << "\t" << *this << std::endl;
#endif
    return *this;
  }

  DilutionIterator operator++() {
    ++block_sink_;

    if (block_sink_ == num_block_) {
      ++block_source_;
      block_sink_ = block_source_;
    }

    return *this;
  }

  bool operator!=(DilutionIterator const &other) const {
    return block_source_ != other.block_source_ || block_sink_ != other.block_sink_ ||
           num_slice_ != other.num_slice_ || num_block_ != other.num_block_;
  }

  BlockIterator begin() const {
    int const block_size = num_slice_ / num_block_;

    if (type_ == DilutionType::block) {
      return BlockIterator(block_source_ * block_size,
                           block_sink_ * block_size,
                           block_source_,
                           block_sink_,
                           num_slice_,
                           num_block_,
                           0,
                           type_,
                           one_sink_slice_);
    } else {
      return BlockIterator(block_source_,
                           block_sink_,
                           block_source_,
                           block_sink_,
                           num_slice_,
                           num_block_,
                           0,
                           type_,
                           one_sink_slice_);
    }
  }

  BlockIterator end() const {
    int const block_size = num_slice_ / num_block_;

    if (type_ == DilutionType::block) {
      return BlockIterator(block_source_ * block_size,
                           block_sink_ * block_size,
                           block_source_,
                           block_sink_,
                           num_slice_,
                           num_block_,
                           2,
                           type_,
                           one_sink_slice_);
    } else {
      return BlockIterator(block_source_,
                           block_sink_,
                           block_source_,
                           block_sink_,
                           num_slice_,
                           num_block_,
                           2,
                           type_,
                           one_sink_slice_);
    }
  }

  int source() const { return block_source_; }

  int sink() const { return block_sink_; }

  DilutionIterator one_sink_slice() const {
    DilutionIterator copy = *this;
    copy.one_sink_slice_ = true;
    return copy;
  }

 private:
  int block_source_;
  int block_sink_;
  int num_slice_;
  int num_block_;
  DilutionType type_;
  bool one_sink_slice_ = false;
};

class DilutionScheme {
 public:
  static DilutionScheme make_full_dilution(int const num_slice) {
    return DilutionScheme{num_slice, num_slice, DilutionType::block};
  }

  DilutionScheme(int const num_slice, int const block_size, DilutionType const type)
      : num_slice_(num_slice), num_block_(num_slice / block_size), type_(type) {}

  DilutionIterator operator[](int const i) const {
    int block_sink = 0;
    int block_source = 0;
    for (int j = 0; j < i; ++j) {
      ++block_sink;

      if (block_sink == num_block_) {
        ++block_source;
        block_sink = block_source;
      }
    }

    return DilutionIterator(block_source, block_sink, num_slice_, num_block_, type_);
  }

  int size() const { return num_block_ * (num_block_ + 1) / 2; }

  DilutionIterator begin() const {
    return DilutionIterator(0, 0, num_slice_, num_block_, type_);
  }

  DilutionIterator end() const {
    return DilutionIterator(num_block_, num_block_, num_slice_, num_block_, type_);
  }

  int time_to_block(int const time) const {
    auto const block_size = num_slice_ / num_block_;

    if (type_ == DilutionType::block) {
      return time / block_size;
    } else if (type_ == DilutionType::interlace) {
      return time % block_size;
    } else {
      throw std::domain_error("This dilution scheme is not implemented.");
    }
  }

 private:
  int num_slice_;
  int num_block_;
  DilutionType type_;
};

void test_dilution_scheme(int const num_slice,
                          int const num_block,
                          DilutionType const type);
