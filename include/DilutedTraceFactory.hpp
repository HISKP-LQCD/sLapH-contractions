#pragma once

#include "DilutedFactorY.hpp"

#include <boost/multi_array.hpp>

template <DilutedFactorType qlt, size_t rvecs>
struct DilutedTrace1Factory {
  static constexpr int num_times = DilutedFactorTypeTraits<qlt>::num_times - 1;

  using Key = std::array<int, num_times>;
  using Value = DilutedTraces<rvecs>;

  DilutedTrace1Factory(DilutedFactorFactory<qlt> &_df,
                      std::vector<Indices> const &_dic,
                      DilutionScheme const &_ds)
      : df(_df), diagram_index_collection(_dic), dilution_scheme(_ds) {}

  Value const &operator[](Key const &key) {
    if (Tr.count(key) == 0) {
      build(key);
    }

    return Tr.at(key);
  }

  void build(Key const &time_key);

  void clear() { return; }

 private:
  DilutedFactorFactory<qlt> &df;
  std::vector<Indices> const &diagram_index_collection;
  DilutionScheme const &dilution_scheme;
  std::map<Key, Value> Tr;
};

template <DilutedFactorType qlt1, DilutedFactorType qlt2, size_t rvecs>
struct DilutedTrace2Factory {
 public:
  /** num_times is the sum of times of contained factors -1 for each continuity
   *  condition of the quarkline diagram
   */
  static constexpr int num_times = DilutedFactorTypeTraits<qlt1>::num_times +
                                   DilutedFactorTypeTraits<qlt2>::num_times - 2;

  using Key = std::array<int, num_times>;
  using Value = DilutedTraces<rvecs>;

  DilutedTrace2Factory(DilutedFactorFactory<qlt1> &_df1,
                           DilutedFactorFactory<qlt2> &_df2,
                           std::vector<Indices> const &_dic,
                           DilutionScheme const &_ds)
      : df1(_df1), df2(_df2), diagram_index_collection(_dic), dilution_scheme(_ds) {}

  Value const &operator[](Key const &key) {
    if (Tr.count(key) == 0) {
      build(key);
    }

    return Tr.at(key);
  }

  void build(Key const &time_key);

  void clear() { Tr.clear(); }

 private:
  DilutedFactorFactory<qlt1> &df1;
  DilutedFactorFactory<qlt2> &df2;
  std::vector<Indices> const &diagram_index_collection;
  DilutionScheme const &dilution_scheme;
  std::map<Key, Value> Tr;
};

template <DilutedFactorType qlt1,
          DilutedFactorType qlt2,
          DilutedFactorType qlt3,
          size_t rvecs>
struct DilutedTrace3Factory {
 public:
  /** num_times is the sum of times of contained factors -1 for each continuity
   *  condition of the quarkline diagram
   */
  static constexpr int num_times = DilutedFactorTypeTraits<qlt1>::num_times +
                                   DilutedFactorTypeTraits<qlt2>::num_times +
                                   DilutedFactorTypeTraits<qlt3>::num_times - 2;

  using Key = std::array<int, num_times>;
  using Value = DilutedTraces<rvecs>;

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
    if (Tr.count(key) == 0) {
      build(key);
    }

    return Tr.at(key);
  }

  void build(Key const &time_key);

  void clear() { Tr.clear(); }

 private:
  DilutedFactorFactory<qlt1> &df1;
  DilutedFactorFactory<qlt2> &df2;
  DilutedFactorFactory<qlt3> &df3;
  std::vector<Indices> const &diagram_index_collection;
  DilutionScheme const &dilution_scheme;
  std::map<Key, Value> Tr;
};
