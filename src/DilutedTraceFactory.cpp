#include "DilutedTraceFactory.hpp"

#include "timings.hpp"

int get_time(BlockIterator const &slice_pair, Location const loc) {
  if (loc == Location::source) {
    return slice_pair.source();
  } else {
    return slice_pair.sink();
  }
}

template <>
void DilutedTrace1Factory<DilutedFactorType::Q1>::build(Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace1Factory<Q1>::build");

  auto t = time_key[0];
  auto b = dilution_scheme.time_to_block(t);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t}][i] = factor_to_trace(df[{t, b}].at({c_look[0]}));
  }
}

template class DilutedTrace1Factory<DilutedFactorType::Q1>;

template <>
void DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1>::build(
    Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace2Factory<Q1>::build");

  auto t1 = time_key[0];
  auto t2 = time_key[1];
  auto b1 = dilution_scheme.time_to_block(t1);
  auto b2 = dilution_scheme.time_to_block(t2);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t1, t2}][i] =
        factor_to_trace(df1[{t1, b2}].at({c_look[0]}), df2[{t2, b1}].at({c_look[1]}));
  }
}

template class DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2>;

template <>
void DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2>::build(
    Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace2Factory<Q0, Q2>::build");

  auto t1 = time_key[0];
  auto t2 = time_key[1];
  auto b2 = dilution_scheme.time_to_block(t2);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t1, t2}][i] =
        factor_to_trace(df1[{t2}].at({c_look[1]}), df2[{b2, t1, b2}].at({c_look[0]}));
  }
}

template class DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1>;

template <>
void DilutedTrace3Factory<DilutedFactorType::Q1,
                          DilutedFactorType::Q1,
                          DilutedFactorType::Q1>::build(Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace3Factory<Q1, Q1, Q1>::build");

  auto const t1 = time_key[0];
  auto const t2 = time_key[1];
  auto const t3 = time_key[2];
  auto const b1 = dilution_scheme.time_to_block(t1);
  auto const b2 = dilution_scheme.time_to_block(t2);
  auto const b3 = dilution_scheme.time_to_block(t3);

  DilutedFactorsMap<2> L1;
  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    multiply<1, 1>(L1, {c_look[0], c_look[1]}, df1[{t1, b2}], df2[{t2, b3}]);
  }

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t1, t2, t3}][i] =
        factor_to_trace(L1[{c_look[0], c_look[1]}], df3[{t3, b1}].at({c_look[2]}));
  }
}

template class DilutedTrace3Factory<DilutedFactorType::Q1,
                                    DilutedFactorType::Q1,
                                    DilutedFactorType::Q1>;

template <>
void DilutedTrace3Factory<DilutedFactorType::Q1,
                          DilutedFactorType::Q0,
                          DilutedFactorType::Q2>::build(Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace3Factory<Q1, Q0, Q2>::build");

  auto const t1 = time_key[0];
  auto const t2 = time_key[1];
  auto const t3 = time_key[2];
  auto const b1 = dilution_scheme.time_to_block(t1);
  auto const b2 = dilution_scheme.time_to_block(t2);
  auto const b3 = dilution_scheme.time_to_block(t3);

  DilutedFactorsMap<2> L1;
  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    multiply<1, 1>(L1, {c_look[2], c_look[0]}, df2[{t1}], df3[{b1, t1, b2}]);
  }

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t1, t2, t3}][i] =
        factor_to_trace(L1[{c_look[2], c_look[0]}], df1[{t2, b3}].at({c_look[1]}));
  }
}

template class DilutedTrace3Factory<DilutedFactorType::Q1,
                                    DilutedFactorType::Q0,
                                    DilutedFactorType::Q2>;

template <>
void DilutedTrace4Factory<DilutedFactorType::Q1,
                          DilutedFactorType::Q1,
                          DilutedFactorType::Q1,
                          DilutedFactorType::Q1>::build(Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace4Factory<Q1, Q1, Q1, Q1>::build");

  auto const t0 = time_key[0];
  auto const t1 = time_key[1];
  auto const t2 = time_key[2];
  auto const t3 = time_key[3];
  auto const b0 = dilution_scheme.time_to_block(t0);
  auto const b1 = dilution_scheme.time_to_block(t1);
  auto const b2 = dilution_scheme.time_to_block(t2);
  auto const b3 = dilution_scheme.time_to_block(t3);

  DilutedFactorsMap<2> L1;
  DilutedFactorsMap<2> L2;
  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    multiply<1, 1>(L1, {c_look[0], c_look[1]}, df1[{t0, b1}], df2[{t1, b2}]);
    multiply<1, 1>(L2, {c_look[2], c_look[3]}, df3[{t2, b3}], df4[{t3, b0}]);
  }

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t0, t1, t2, t3}][i] =
        factor_to_trace(L1[{c_look[0], c_look[1]}], L2[{c_look[2], c_look[3]}]);
  }
}

template class DilutedTrace4Factory<DilutedFactorType::Q1,
                                    DilutedFactorType::Q1,
                                    DilutedFactorType::Q1,
                                    DilutedFactorType::Q1>;

template <>
void DilutedTrace4Factory<DilutedFactorType::Q2,
                          DilutedFactorType::Q0,
                          DilutedFactorType::Q2,
                          DilutedFactorType::Q0>::build(Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace4Factory<Q2, Q0, Q2, Q0>::build");

  auto const t0 = time_key[0];
  auto const t1 = time_key[1];
  auto const t2 = time_key[2];
  auto const t3 = time_key[3];
  auto const b0 = dilution_scheme.time_to_block(t0);
  // auto const b1 = dilution_scheme.time_to_block(t1);
  auto const b2 = dilution_scheme.time_to_block(t2);
  // auto const b3 = dilution_scheme.time_to_block(t3);

  DilutedFactorsMap<2> L1;
  DilutedFactorsMap<2> L2;
  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    multiply<1, 1>(L1, {c_look[0], c_look[1]}, df2[{t0}], df3[{b0, t1, b2}]);
    multiply<1, 1>(L2, {c_look[2], c_look[3]}, df4[{t2}], df1[{b2, t3, b0}]);
  }

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[time_key][i] =
        factor_to_trace(L1[{c_look[0], c_look[1]}], L2[{c_look[2], c_look[3]}]);
  }
}

template class DilutedTrace4Factory<DilutedFactorType::Q2,
                                    DilutedFactorType::Q0,
                                    DilutedFactorType::Q2,
                                    DilutedFactorType::Q0>;

template <>
void DilutedTrace6Factory<DilutedFactorType::Q2,
                          DilutedFactorType::Q0,
                          DilutedFactorType::Q2,
                          DilutedFactorType::Q0,
                          DilutedFactorType::Q2,
                          DilutedFactorType::Q0>::build(Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace6Factory<Q2, Q0, Q2, Q0, Q2, Q0>::build");

  auto const t0 = time_key[0];
  auto const t1 = time_key[1];
  auto const t2 = time_key[2];
  auto const t3 = time_key[3];
  auto const t4 = time_key[4];
  auto const t5 = time_key[5];
  // auto const b0 = dilution_scheme.time_to_block(t0);
  auto const b1 = dilution_scheme.time_to_block(t1);
  // auto const b2 = dilution_scheme.time_to_block(t2);
  auto const b3 = dilution_scheme.time_to_block(t3);
  // auto const b4 = dilution_scheme.time_to_block(t4);
  auto const b5 = dilution_scheme.time_to_block(t5);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    auto const &c_look = diagram_index_collection[i];

    auto const &m0 = df1[{b5, t0, b1}];
    auto const &m1 = df2[{t1}];
    auto const &m2 = df3[{b1, t2, b3}];
    auto const &m3 = df4[{t3}];
    auto const &m4 = df5[{b3, t4, b5}];
    auto const &m5 = df6[{t5}];

    auto const &f0 = m0.at({c_look[1]});
    auto const &f1 = m1.at({c_look[2]});
    auto const &f2 = m2.at({c_look[3]});
    auto const &f3 = m3.at({c_look[4]});
    auto const &f4 = m4.at({c_look[5]});
    auto const &f5 = m5.at({c_look[0]});

    Tr[time_key][i] = factor_to_trace(f0 * f1 * f2 * f3 * f4, f5);
  }
}

template class DilutedTrace6Factory<DilutedFactorType::Q2,
                                    DilutedFactorType::Q0,
                                    DilutedFactorType::Q2,
                                    DilutedFactorType::Q0,
                                    DilutedFactorType::Q2,
                                    DilutedFactorType::Q0>;
