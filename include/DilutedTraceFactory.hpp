#pragma once

#include "DilutedFactorFactory.hpp"
#include "DilutedProductFactory.hpp"
#include "DilutedTrace.hpp"

#include <boost/multi_array.hpp>

int get_time(BlockIterator const &slice_pair, Location const loc);

template <int num_times>
std::array<int, num_times> make_key(BlockIterator const &slice_pair,
                                    std::vector<Location> const &locations) {
  std::array<int, num_times> key;
  std::transform(std::begin(locations),
                 std::end(locations),
                 std::begin(key),
                 [&slice_pair](Location const loc) { return get_time(slice_pair, loc); });
  return key;
}

class AbstractDilutedTraceFactory {
 public:
  virtual ~AbstractDilutedTraceFactory() {}

  virtual DilutedTracesMap const &get(BlockIterator const &slice_pair,
                                      std::vector<Location> const &locations) = 0;

  virtual void clear() = 0;
};

template <DilutedFactorType qlt>
class DilutedTrace1Factory : public AbstractDilutedTraceFactory {
 public:
  static constexpr int num_times = DilutedFactorTypeTraits<qlt>::num_times - 1;

  using Key = std::array<int, num_times>;
  using Value = DilutedTracesMap;

  DilutedTrace1Factory(DilutedFactorFactory<qlt> &_df,
                       std::vector<Indices> const &_dic,
                       DilutionScheme const &_ds)
      : df(_df), diagram_index_collection(_dic), dilution_scheme(_ds) {}

  Value const &operator[](Key const &key) {
    Value *result = nullptr;
#pragma omp critical(DilutedTrace1Factory_operator_square)
    {
      if (Tr.count(key) == 0) {
        build(key);
      }

      result = &Tr.at(key);
    }
    return *result;
  }

  DilutedTracesMap const &get(BlockIterator const &slice_pair,
                              std::vector<Location> const &locations) override {
    assert(ssize(locations) == num_times);
    auto const &key = make_key<num_times>(slice_pair, locations);
    return (*this)[key];
  }

  void build(Key const &time_key);

  void clear() override { return; }

 private:
  DilutedFactorFactory<qlt> &df;
  std::vector<Indices> const &diagram_index_collection;
  DilutionScheme const &dilution_scheme;
  std::map<Key, Value> Tr;
};

template <DilutedFactorType qlt1, DilutedFactorType qlt2>
class DilutedTrace2Factory : public AbstractDilutedTraceFactory {
 public:
  /** num_times is the sum of times of contained factors -1 for each continuity
   *  condition of the quarkline diagram
   */
  static constexpr int num_times = DilutedFactorTypeTraits<qlt1>::num_times +
                                   DilutedFactorTypeTraits<qlt2>::num_times - 2;

  using Key = std::array<int, num_times>;
  using Value = DilutedTracesMap;

  DilutedTrace2Factory(DilutedFactorFactory<qlt1> &_df1,
                       DilutedFactorFactory<qlt2> &_df2,
                       std::vector<Indices> const &_dic,
                       DilutionScheme const &_ds)
      : df1(_df1), df2(_df2), diagram_index_collection(_dic), dilution_scheme(_ds) {}

  Value const &operator[](Key const &key) {
    Value *result = nullptr;
#pragma omp critical(DilutedTrace2Factory_operator_square)
    {
      if (Tr.count(key) == 0) {
        build(key);
      }

      result = &Tr.at(key);
    }
    return *result;
  }

  DilutedTracesMap const &get(BlockIterator const &slice_pair,
                              std::vector<Location> const &locations) override {
    assert(ssize(locations) == num_times);
    auto const &key = make_key<num_times>(slice_pair, locations);
    return (*this)[key];
  }

  void build(Key const &time_key);

  void clear() override { Tr.clear(); }

 private:
  DilutedFactorFactory<qlt1> &df1;
  DilutedFactorFactory<qlt2> &df2;
  std::vector<Indices> const &diagram_index_collection;
  DilutionScheme const &dilution_scheme;
  std::map<Key, Value> Tr;
};

template <DilutedFactorType qlt1, DilutedFactorType qlt2, DilutedFactorType qlt3>
class DilutedTrace3Factory : public AbstractDilutedTraceFactory {
 public:
  /** num_times is the sum of times of contained factors -1 for each continuity
   *  condition of the quarkline diagram
   */
  static constexpr int num_times = DilutedFactorTypeTraits<qlt1>::num_times +
                                   DilutedFactorTypeTraits<qlt2>::num_times +
                                   DilutedFactorTypeTraits<qlt3>::num_times - 3;

  using Key = std::array<int, num_times>;
  using Value = DilutedTracesMap;

  DilutedTrace3Factory(DilutedFactorFactory<qlt1> &_df1,
                       DilutedFactorFactory<qlt2> &_df2,
                       DilutedFactorFactory<qlt3> &_df3,
                       std::vector<Indices> const &_dic,
                       DilutionScheme const &_ds)
      : df1(_df1),
        df2(_df2),
        df3(_df3),
        diagram_index_collection(_dic),
        dilution_scheme(_ds) {}

  Value const &operator[](Key const &key) {
    Value *result = nullptr;
#pragma omp critical(DilutedTrace3Factory_operator_square)
    {
      if (Tr.count(key) == 0) {
        build(key);
      }

      result = &Tr.at(key);
    }
    return *result;
  }

  DilutedTracesMap const &get(BlockIterator const &slice_pair,
                              std::vector<Location> const &locations) override {
    assert(ssize(locations) == num_times);
    auto const &key = make_key<num_times>(slice_pair, locations);
    return (*this)[key];
  }

  void build(Key const &time_key);

  void clear() override { Tr.clear(); }

 private:
  DilutedFactorFactory<qlt1> &df1;
  DilutedFactorFactory<qlt2> &df2;
  DilutedFactorFactory<qlt3> &df3;
  std::vector<Indices> const &diagram_index_collection;
  DilutionScheme const &dilution_scheme;
  std::map<Key, Value> Tr;
};

template <DilutedFactorType qlt1,
          DilutedFactorType qlt2,
          DilutedFactorType qlt3,
          DilutedFactorType qlt4>
class DilutedTrace4Factory : public AbstractDilutedTraceFactory {
 public:
  /** num_times is the sum of times of contained factors -1 for each continuity
   *  condition of the quarkline diagram
   */
  static constexpr int num_times = DilutedFactorTypeTraits<qlt1>::num_times +
                                   DilutedFactorTypeTraits<qlt2>::num_times +
                                   DilutedFactorTypeTraits<qlt3>::num_times +
                                   DilutedFactorTypeTraits<qlt4>::num_times - 4;

  using Key = std::array<int, num_times>;
  using Value = DilutedTracesMap;

  DilutedTrace4Factory(DilutedFactorFactory<qlt1> &_df1,
                       DilutedFactorFactory<qlt2> &_df2,
                       DilutedFactorFactory<qlt3> &_df3,
                       DilutedFactorFactory<qlt4> &_df4,
                       DilutedProductFactoryQ0Q2 &dpf,
                       std::vector<Indices> const &_dic,
                       DilutionScheme const &_ds)
      : df1(_df1),
        df2(_df2),
        df3(_df3),
        df4(_df4),
        dpf_(dpf),
        diagram_index_collection(_dic),
        dilution_scheme(_ds) {}

  Value const &operator[](Key const &key) {
    Value *result = nullptr;
#pragma omp critical(DilutedTrace4Factory_operator_square)
    {
      if (Tr.count(key) == 0) {
        build(key);
      }

      result = &Tr.at(key);
    }
    return *result;
  }

  DilutedTracesMap const &get(BlockIterator const &slice_pair,
                              std::vector<Location> const &locations) override {
    assert(ssize(locations) == num_times);
    auto const &key = make_key<num_times>(slice_pair, locations);
    return (*this)[key];
  }

  void build(Key const &time_key);

  void clear() override { Tr.clear(); }

 private:
  DilutedFactorFactory<qlt1> &df1;
  DilutedFactorFactory<qlt2> &df2;
  DilutedFactorFactory<qlt3> &df3;
  DilutedFactorFactory<qlt4> &df4;
  DilutedProductFactoryQ0Q2 &dpf_;
  std::vector<Indices> const &diagram_index_collection;
  DilutionScheme const &dilution_scheme;
  std::map<Key, Value> Tr;
};

template <DilutedFactorType qlt1,
          DilutedFactorType qlt2,
          DilutedFactorType qlt3,
          DilutedFactorType qlt4,
          DilutedFactorType qlt5,
          DilutedFactorType qlt6>
class DilutedTrace6Factory : public AbstractDilutedTraceFactory {
 public:
  /** num_times is the sum of times of contained factors -1 for each continuity
   *  condition of the quarkline diagram
   */
  static constexpr int num_times = DilutedFactorTypeTraits<qlt1>::num_times +
                                   DilutedFactorTypeTraits<qlt2>::num_times +
                                   DilutedFactorTypeTraits<qlt3>::num_times +
                                   DilutedFactorTypeTraits<qlt4>::num_times +
                                   DilutedFactorTypeTraits<qlt5>::num_times +
                                   DilutedFactorTypeTraits<qlt6>::num_times - 6;

  using Key = std::array<int, num_times>;
  using Value = DilutedTracesMap;

  DilutedTrace6Factory(DilutedFactorFactory<qlt1> &_df1,
                       DilutedFactorFactory<qlt2> &_df2,
                       DilutedFactorFactory<qlt3> &_df3,
                       DilutedFactorFactory<qlt4> &_df4,
                       DilutedFactorFactory<qlt5> &_df5,
                       DilutedFactorFactory<qlt6> &_df6,
                       DilutedProductFactoryQ0Q2 &dpf,
                       std::vector<Indices> const &_dic,
                       DilutionScheme const &_ds)
      : df1(_df1),
        df2(_df2),
        df3(_df3),
        df4(_df4),
        df5(_df5),
        df6(_df6),
        dpf_(dpf),
        diagram_index_collection(_dic),
        dilution_scheme(_ds) {}

  Value const &operator[](Key const &key) {
    Value *result = nullptr;
#pragma omp critical(DilutedTrace6Factory_operator_square)
    {
      if (Tr.count(key) == 0) {
        build(key);
      }

      result = &Tr.at(key);
    }
    return *result;
  }

  DilutedTracesMap const &get(BlockIterator const &slice_pair,
                              std::vector<Location> const &locations) override {
    assert(ssize(locations) == num_times);
    auto const &key = make_key<num_times>(slice_pair, locations);
    return (*this)[key];
  }

  void build(Key const &time_key);

  void clear() override { Tr.clear(); }

 private:
  DilutedFactorFactory<qlt1> &df1;
  DilutedFactorFactory<qlt2> &df2;
  DilutedFactorFactory<qlt3> &df3;
  DilutedFactorFactory<qlt4> &df4;
  DilutedFactorFactory<qlt5> &df5;
  DilutedFactorFactory<qlt6> &df6;
  DilutedProductFactoryQ0Q2 &dpf_;
  std::vector<Indices> const &diagram_index_collection;
  DilutionScheme const &dilution_scheme;
  std::map<Key, Value> Tr;
};
