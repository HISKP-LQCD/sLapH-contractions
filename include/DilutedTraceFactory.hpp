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

template <DilutedFactorType qlt, DilutedFactorType... qlts>
constexpr int get_num_times() {
  return DilutedFactorTypeTraits<qlt>::num_times - 1 + get_num_times<qlts...>;
}

template <DilutedFactorType qlt>
constexpr int get_num_times() {
  return DilutedFactorTypeTraits<qlt>::num_times - 1;
}

class AbstractDilutedTraceFactory {
 public:
  virtual ~AbstractDilutedTraceFactory() {}

  virtual void request(BlockIterator const &slice_pair,
                       std::vector<Location> const &locations) {
    abort();
  }

  virtual void build_all() { abort(); }

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
    //#pragma omp critical(DilutedTrace1Factory_operator_square)
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
    //#pragma omp critical(DilutedTrace2Factory_operator_square)
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
    //#pragma omp critical(DilutedTrace3Factory_operator_square)
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

  void request(Key const &time_key) {
    // We store the request for the trace.
    requests_.insert(time_key);

    // And we can also schedule the `DilutedFactoryFactory` from here.
    auto const t0 = time_key[0];
    auto const t1 = time_key[1];
    auto const t2 = time_key[2];
    auto const t3 = time_key[3];
    auto const b0 = dilution_scheme.time_to_block(t0);
    auto const b2 = dilution_scheme.time_to_block(t2);

    for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
      const auto &c_look = diagram_index_collection[i];
      dpf_.request({b0, t1, b2, t2}, {c_look[1], c_look[2]});
      dpf_.request({b2, t3, b0, t0}, {c_look[3], c_look[0]});
    }
  }

  void build_all() {
    // We call `build_all` on the `DilutedProductFactory` here, but we want to pull that
    // out further such that we can request more before starting to build.
    dpf_.build_all();

    std::vector<Key> unique_requests;
    unique_requests.reserve(requests_.size());
    for (auto const &time_key : requests_) {
      if (Tr.count(time_key) == 0) {
        unique_requests.push_back(time_key);
        Tr[time_key];
      }
    }

    requests_.clear();

    for (auto i = 0; i < ssize(unique_requests); ++i) {
      build(unique_requests[i]);
    }
  }

  Value const &operator[](Key const &key) { return Tr.at(key); }

  DilutedTracesMap const &get(BlockIterator const &slice_pair,
                              std::vector<Location> const &locations) override {
    assert(ssize(locations) == num_times);
    auto const &key = make_key<num_times>(slice_pair, locations);
    return (*this)[key];
  }

  void clear() override { Tr.clear(); }

 private:
  void build(Key const &time_key);

  std::set<Key> requests_;

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
    //#pragma omp critical(DilutedTrace6Factory_operator_square)
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
