#include "DilutedTraceFactory.hpp"

#include "timings.hpp"

int get_time(BlockIterator const &slice_pair, Location const loc) {
  if (loc == Location::source) {
    return slice_pair.source();
  } else {
    return slice_pair.sink();
  }
}

/*****************************************************************************/
/*                                    Q1                                     */
/*****************************************************************************/

template <>
void DilutedTrace1Factory<DilutedFactorType::Q1>::build(Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace1Factory<Q1>::build");

  auto t = time_key[0];
  auto b = dilution_scheme.time_to_block(t);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    auto const &value = factor_to_trace(df[{t, b}].at({c_look[0]}));
    Tr[{t}][i] = value;
  }
}

template class DilutedTrace1Factory<DilutedFactorType::Q1>;

/*****************************************************************************/
/*                                   Q1 Q1                                   */
/*****************************************************************************/

template <>
void DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1>::build(
    Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace2Factory<Q1>::build");

  auto t1 = time_key[0];
  auto t2 = time_key[1];
  auto b1 = dilution_scheme.time_to_block(t1);
  auto b2 = dilution_scheme.time_to_block(t2);

  df1.request({t1, b2});
  df2.request({t2, b1});
  df1.build_all();
  df2.build_all();

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    auto const &value =
        factor_to_trace(df1[{t1, b2}].at({c_look[0]}), df2[{t2, b1}].at({c_look[1]}));
    Tr[{t1, t2}][i] = value;
  }
}

template class DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2>;

/*****************************************************************************/
/*                                   Q0 Q2                                   */
/*****************************************************************************/

template <>
void DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2>::build(
    Key const &time_key) {
  TimingScope<3> timing_scope("DilutedTrace2Factory<Q0, Q2>::build");

  auto t1 = time_key[0];
  auto t2 = time_key[1];
  auto b2 = dilution_scheme.time_to_block(t2);

  df1.request({t2});
  df2.request({b2, t1, b2});
  df1.build_all();
  df2.build_all();

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    auto const &c_look = diagram_index_collection[i];
    auto const &value =
        factor_to_trace(df1[{t2}].at({c_look[1]}), df2[{b2, t1, b2}].at({c_look[0]}));
    Tr[{t1, t2}][i] = value;
  }
}

template class DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1>;

/*****************************************************************************/
/*                                 Q1 Q1 Q1                                  */
/*****************************************************************************/

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

  df1.request({t1, b2});
  df2.request({t2, b3});
  df3.request({t3, b1});
  df1.build_all();
  df2.build_all();
  df3.build_all();

  DilutedFactorsMap<2> L1;
  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    multiply<1, 1>(L1, {c_look[0], c_look[1]}, df1[{t1, b2}], df2[{t2, b3}]);
  }

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    auto const &value =
        factor_to_trace(L1[{c_look[0], c_look[1]}], df3[{t3, b1}].at({c_look[2]}));
    Tr[{t1, t2, t3}][i] = value;
  }
}

template class DilutedTrace3Factory<DilutedFactorType::Q1,
                                    DilutedFactorType::Q1,
                                    DilutedFactorType::Q1>;

/*****************************************************************************/
/*                                 Q1 Q0 Q2                                  */
/*****************************************************************************/

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

  df1.request({t2, b3});
  df2.request({t1});
  df3.request({b1, t1, b2});
  df1.build_all();
  df2.build_all();
  df3.build_all();

  DilutedFactorsMap<2> L1;
  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    multiply<1, 1>(L1, {c_look[2], c_look[0]}, df2[{t1}], df3[{b1, t1, b2}]);
  }

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    auto const &value =
        factor_to_trace(L1[{c_look[2], c_look[0]}], df1[{t2, b3}].at({c_look[1]}));
    Tr[{t1, t2, t3}][i] = value;
  }
}

template class DilutedTrace3Factory<DilutedFactorType::Q1,
                                    DilutedFactorType::Q0,
                                    DilutedFactorType::Q2>;

/*****************************************************************************/
/*                                Q1 Q1 Q1 Q1                                */
/*****************************************************************************/

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

  df1.request({t0, b1});
  df2.request({t1, b2});
  df3.request({t2, b3});
  df4.request({t3, b0});
  df1.build_all();
  df2.build_all();
  df3.build_all();
  df4.build_all();

  DilutedFactorsMap<2> L1;
  DilutedFactorsMap<2> L2;
  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    multiply<1, 1>(L1, {c_look[0], c_look[1]}, df1[{t0, b1}], df2[{t1, b2}]);
    multiply<1, 1>(L2, {c_look[2], c_look[3]}, df3[{t2, b3}], df4[{t3, b0}]);
  }

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    auto const &value =
        factor_to_trace(L1[{c_look[0], c_look[1]}], L2[{c_look[2], c_look[3]}]);
    Tr[{t0, t1, t2, t3}][i] = value;
  }
}

template class DilutedTrace4Factory<DilutedFactorType::Q1,
                                    DilutedFactorType::Q1,
                                    DilutedFactorType::Q1,
                                    DilutedFactorType::Q1>;

/*****************************************************************************/
/*                                Q2 Q0 Q2 Q0                                */
/*****************************************************************************/

template <>
void DilutedTrace4Factory<DilutedFactorType::Q2,
                          DilutedFactorType::Q0,
                          DilutedFactorType::Q2,
                          DilutedFactorType::Q0>::request_impl(Key const &time_key) {
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
  auto const b2 = dilution_scheme.time_to_block(t2);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];

    auto const &l1 = dpf_.get({b0, t1, b2, t2}, {c_look[1], c_look[2]});
    auto const &l2 = dpf_.get({b2, t3, b0, t0}, {c_look[3], c_look[0]});
    auto const &value = factor_to_trace(l1, l2);

    Tr[time_key][i] = value;
  }

#ifdef SLAPH_CLEAR_QQ_CACHE
  dpf_.clear();
#endif
}

template class DilutedTrace4Factory<DilutedFactorType::Q2,
                                    DilutedFactorType::Q0,
                                    DilutedFactorType::Q2,
                                    DilutedFactorType::Q0>;

/*****************************************************************************/
/*                             Q2 Q0 Q2 Q0 Q2 Q0                             */
/*****************************************************************************/

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
  auto const b0 = dilution_scheme.time_to_block(t0);
  auto const b2 = dilution_scheme.time_to_block(t2);
  auto const b4 = dilution_scheme.time_to_block(t4);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    auto const &c_look = diagram_index_collection[i];

    dpf_.request({b0, t1, b2, t2}, {c_look[1], c_look[2]});
    dpf_.request({b2, t3, b4, t4}, {c_look[3], c_look[4]});
    dpf_.request({b4, t5, b0, t0}, {c_look[5], c_look[0]});
  }

  dpf_.build_all();

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    auto const &c_look = diagram_index_collection[i];

    auto const &l01 = dpf_.get({b0, t1, b2, t2}, {c_look[1], c_look[2]});
    auto const &l23 = dpf_.get({b2, t3, b4, t4}, {c_look[3], c_look[4]});
    auto const &l45 = dpf_.get({b4, t5, b0, t0}, {c_look[5], c_look[0]});

    assert(l01.size() > 0);
    assert(l23.size() > 0);
    assert(l45.size() > 0);
    auto const &value = factor_to_trace(l01 * l23, l45);
    Tr[time_key][i] = value;
  }

#ifdef SLAPH_CLEAR_QQ_CACHE
  dpf_.clear();
#endif
}

template class DilutedTrace6Factory<DilutedFactorType::Q2,
                                    DilutedFactorType::Q0,
                                    DilutedFactorType::Q2,
                                    DilutedFactorType::Q0,
                                    DilutedFactorType::Q2,
                                    DilutedFactorType::Q0>;
