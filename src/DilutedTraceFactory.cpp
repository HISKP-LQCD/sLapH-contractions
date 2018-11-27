#include "DilutedTraceFactory.hpp"

template <>
void DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>::build(
    Key const &time_key) {
  auto t1 = time_key[0];
  auto t2 = time_key[1];
  auto b2 = dilution_scheme.time_to_block(t2);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t1, t2}][i] = factor_to_trace(df1[{t2}].at({c_look[1]}),
                                      df2[{b2, t1, b2}].at({c_look[0]}));
  }
}

template <>
void DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>::build(
    Key const &time_key) {
  auto t1 = time_key[0];
  auto t2 = time_key[1];
  auto b1 = dilution_scheme.time_to_block(t1);
  auto b2 = dilution_scheme.time_to_block(t2);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t1, t2}][i] = factor_to_trace(df1[{t1, b2}].at({c_look[0]}),
                                      df2[{t2, b1}].at({c_look[1]}));
  }
}

template <>
void DilutedTrace1Factory<DilutedFactorType::Q1, 1>::build(Key const &time_key) {
  auto t = time_key[0];
  auto b = dilution_scheme.time_to_block(t);

  for (ssize_t i = 0; i != ssize(diagram_index_collection); ++i) {
    const auto &c_look = diagram_index_collection[i];
    Tr[{t}][i] = factor_to_trace(df[{t, b}].at({c_look[0]}));
  }
}

template class DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>;
template class DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>;
template class DilutedTrace1Factory<DilutedFactorType::Q1, 1>;
