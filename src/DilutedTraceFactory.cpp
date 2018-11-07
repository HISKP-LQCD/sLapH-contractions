#include "DilutedTraceFactory.h"

template <>
void DilutedTraceTraceFactory<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>::build(
    Key const &time_key) {
  auto t1 = time_key[0];
  auto t2 = time_key[1];
  auto b2 = dilution_scheme.time_to_block(t2);

  for (const auto &c_look : diagram_index_collection) {
    Tr[{t1, t2}][c_look.id] = factor_to_trace(df1[{t2}].at({c_look.lookup[1]}),
                                              df2[{b2, t1, b2}].at({c_look.lookup[0]}));
  }
}

template <>
void DilutedTraceTraceFactory<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>::build(
    Key const &time_key) {
  auto t1 = time_key[0];
  auto t2 = time_key[1];
  auto b1 = dilution_scheme.time_to_block(t1);
  auto b2 = dilution_scheme.time_to_block(t2);

  for (const auto &c_look : diagram_index_collection) {
    Tr[{t1, t2}][c_look.id] = factor_to_trace(df1[{t1, b2}].at({c_look.lookup[0]}),
                                              df2[{t2, b1}].at({c_look.lookup[1]}));
  }
}

template class DilutedTraceTraceFactory<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>;
template class DilutedTraceTraceFactory<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>;

template <>
void DilutedTraceFactory<DilutedFactorType::Q1, 1>::build(Key const &time_key) {
  auto t = time_key[0];
  auto b = dilution_scheme.time_to_block(t);

  for (const auto &c_look : diagram_index_collection) {
    Tr[{t}][c_look.id] = factor_to_trace(df[{t, b}].at({c_look.lookup[0]}));
  }
}

template class DilutedTraceFactory<DilutedFactorType::Q1, 1>;
