#include "Diagram.h"

#include <omp.h>

#include "local_timer.h"

/*****************************************************************************/
/*                                    C2c                                    */
/*****************************************************************************/

void C2c::assemble_impl(std::vector<Complex> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    auto const &x = q.trQ0Q2[{slice_pair.source(),slice_pair.sink()}].at(c_look.lookup[0]);
    c[i] += std::accumulate(std::begin(x), std::end(x), Complex(0.0, 0.0)) /
            static_cast<double>(x.size());
  }
}

/*****************************************************************************/
/*                                    C20                                    */
/*****************************************************************************/

void C20::assemble_impl(std::vector<Complex> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    auto const &x = q.trQ1Q1[{slice_pair.source(),slice_pair.sink()}].at(c_look.lookup[0]);
    c[i] += std::accumulate(std::begin(x), std::end(x), Complex(0.0, 0.0)) /
            static_cast<double>(x.size());
  }
}

/*****************************************************************************/
/*                                    C20V                                   */
/*****************************************************************************/

void C20V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] += inner_product(q.trQ1[{slice_pair.source()}].at(c_look.lookup[0]), 
                          q.trQ1[{slice_pair.sink()}].at(c_look.lookup[1]));
  }
}

/*****************************************************************************/
/*                                    C3c                                    */
/*****************************************************************************/

C3c::C3c(std::vector<DiagramIndex> const &corr_lookup,
         std::string const &output_path,
         std::string const &output_filename,
         int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        make_tuple(std::array<ssize_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<ssize_t, 1>{c_look.lookup[1]}));
  }
}

void C3c::assemble_impl(std::vector<Complex> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  DilutedFactors<2, 1> L1;

  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(
        L1,
        std::get<0>(ids),
        q.q0[{slice_pair.source()}],
        q.q2l[{slice_pair.source_block(), slice_pair.source(), slice_pair.sink_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C3c::assemble_impl] multiply");
  
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] +=
        trace(L1[std::get<0>(ids)],
              q.q1[{slice_pair.sink(), slice_pair.source_block()}].at(std::get<1>(ids)));
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C3c::assemble_impl] trace");
}

/*****************************************************************************/
/*                                    C30                                    */
/*****************************************************************************/

C30::C30(std::vector<DiagramIndex> const &corr_lookup,
         std::string const &output_path,
         std::string const &output_filename,
         int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        make_tuple(std::array<ssize_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<ssize_t, 1>{c_look.lookup[1]}));
  }
}

void C30::assemble_impl(std::vector<Complex> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  DilutedFactors<2, 1> L1;
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(L1,
                         std::get<0>(ids),
                         q.q1[{slice_pair.source(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C30::assemble_impl] multiply");

  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] +=
        trace(L1[std::get<0>(ids)],
              q.q1[{slice_pair.sink(), slice_pair.source_block()}].at(std::get<1>(ids)));
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C30::assemble_impl] trace");
}

/*****************************************************************************/
/*                                   C30V                                    */
/*****************************************************************************/

void C30V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

//    assert(c_look.lookup[0] < q.trQ1Q1.tr.shape()[0]);
//    assert(slice_pair.source() < q.trQ1Q1.tr.shape()[1]);
//    assert(slice_pair.source() < q.trQ1Q1.tr.shape()[2]);
//
//    assert(c_look.lookup[1] < q.trQ1.tr.shape()[0]);
//    assert(slice_pair.sink() < q.trQ1.tr.shape()[1]);

    c[i] += inner_product(
        q.trQ1Q1[{slice_pair.source(),slice_pair.source()}].at(c_look.lookup[0]),
        q.trQ1[{slice_pair.sink()}].at(c_look.lookup[1]));
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C30::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C4cD                                    */
/*****************************************************************************/

void C4cD::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] +=
        inner_product(q.trQ0Q2[{slice_pair.source(),slice_pair.sink()}].at(c_look.lookup[0]),
                      q.trQ0Q2[{slice_pair.source(),slice_pair.sink()}].at(c_look.lookup[1]));
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cD::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C40D                                    */
/*****************************************************************************/

void C40D::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] +=
        inner_product(q.trQ1Q1[{slice_pair.source(),slice_pair.sink()}].at(c_look.lookup[0]),
                      q.trQ1Q1[{slice_pair.source(),slice_pair.sink()}].at(c_look.lookup[1]));
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40D::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C4cV                                    */
/*****************************************************************************/

void C4cV::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] += inner_product(
        q.trQ0Q2[{slice_pair.source(),slice_pair.source()}].at(c_look.lookup[0]),
        q.trQ0Q2[{slice_pair.sink(),slice_pair.sink()}].at(c_look.lookup[1]));
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cV::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C40V                                    */
/*****************************************************************************/

void C40V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] += inner_product(
        q.trQ1Q1[{slice_pair.source(),slice_pair.source()}].at(c_look.lookup[0]),
        q.trQ1Q1[{slice_pair.sink(),slice_pair.sink()}].at(c_look.lookup[1]));
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40V::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C4cB                                    */
/*****************************************************************************/

C4cB::C4cB(std::vector<DiagramIndex> const &corr_lookup,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C4cB::assemble_impl(std::vector<Complex> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  DilutedFactors<2, 1> L1;
  DilutedFactors<2, 1> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(
        L1,
        ids[0],
        q.q0[{slice_pair.source()}],
        q.q2l[{slice_pair.source_block(), slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1, 0, 0>(
        L2,
        ids[1],
        q.q0[{slice_pair.sink()}],
        q.q2l[{slice_pair.sink_block(), slice_pair.sink(), slice_pair.source_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cB::assemble_impl] multiply");
  
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cB::assemble_impl] trace");
}

/*****************************************************************************/
/*                                   C40B                                    */
/*****************************************************************************/

C40B::C40B(std::vector<DiagramIndex> const &corr_lookup,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C40B::assemble_impl(std::vector<Complex> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  DilutedFactors<2, 1> L1;
  DilutedFactors<2, 1> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(L1,
                         ids[0],
                         q.q1[{slice_pair.source(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1, 0, 0>(L2,
                         ids[1],
                         q.q1[{slice_pair.sink(), slice_pair.sink_block()}],
                         q.q1[{slice_pair.sink(), slice_pair.source_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40B::assemble_impl] multiply");
  
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40B::assemble_impl] trace");
}

/*****************************************************************************/
/*                                   C4cC                                    */
/*****************************************************************************/

C4cC::C4cC(std::vector<DiagramIndex> const &corr_lookup,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C4cC::assemble_impl(std::vector<Complex> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  DilutedFactors<2, 1> L1;
  DilutedFactors<2, 1> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(
        L1,
        ids[0],
        q.q0[{slice_pair.sink()}],
        q.q2v[{slice_pair.sink_block(), slice_pair.source(), slice_pair.sink_block()}]);
    multiply<1, 1, 0, 0>(
        L2,
        ids[1],
        q.q0[{slice_pair.sink()}],
        q.q2v[{slice_pair.sink_block(), slice_pair.source(), slice_pair.sink_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cC::assemble_impl] multiply");
  
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cC::assemble_impl] trace");
}

/*****************************************************************************/
/*                                   C40C                                    */
/*****************************************************************************/

C40C::C40C(std::vector<DiagramIndex> const &corr_lookup,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C40C::assemble_impl(std::vector<Complex> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  DilutedFactors<2, 1> L1;
  DilutedFactors<2, 1> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(L1,
                         ids[0],
                         q.q1[{slice_pair.sink(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1, 0, 0>(L2,
                         ids[1],
                         q.q1[{slice_pair.sink(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40C::assemble_impl] multiply");
  
  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40C::assemble_impl] trace");
}
