#pragma once

#include "Correlators.hpp"
#include "DilutedFactorY.hpp"
#include "DilutedTraceFactory.hpp"
#include "h5-wrapper.hpp"

#include <omp.h>

#include <mutex>

struct DiagramParts {
  DiagramParts(RandomVector const &random_vector,
               Perambulator const &perambulator,
               OperatorFactory const &meson_operator,
               DilutionScheme const &dilution_scheme,
               ssize_t const dilT,
               ssize_t const dilE,
               ssize_t const nev,
               ssize_t const Lt,
               DilutedFactorIndicesCollection const &dil_fac_lookup,
               DiagramIndicesCollection const &corr_lookup,
               TraceIndicesCollection const &trace_indices_map)
      : q0(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q0")),
        q1(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q1")),
        q2(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q2")),
        trQ1(q1, trace_indices_map.at("trQ1"), dilution_scheme),
        trQ1Q1(q1, q1, trace_indices_map.at("trQ1Q1"), dilution_scheme),
        trQ0Q2(q0, q2, trace_indices_map.at("trQ0Q2"), dilution_scheme),
        trQ1Q1Q1(q1, q1, q1, trace_indices_map.at("trQ1Q1Q1"), dilution_scheme) {
    trace_factories["trQ1"] = &trQ1;
    trace_factories["trQ1Q1"] = &trQ1Q1;
    trace_factories["trQ0Q2"] = &trQ0Q2;
    trace_factories["trQ1Q1Q1"] = &trQ1Q1Q1;
  }

  void clear() {
    q0.clear();
    q1.clear();
    q2.clear();

    for (auto const &elem : trace_factories) {
      elem.second->clear();
    }
  }

  DilutedFactorFactory<DilutedFactorType::Q0> q0;
  DilutedFactorFactory<DilutedFactorType::Q1> q1;
  DilutedFactorFactory<DilutedFactorType::Q2> q2;

  //< Temporal memory for tr(Q1)
  DilutedTrace1Factory<DilutedFactorType::Q1> trQ1;

  //< Temporal memory for tr(rVdaggerV*Q1*rVdaggerV*Q1)
  DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1> trQ1Q1;

  //< Temporal memory for tr(Q2*rVdaggerVr)
  DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2> trQ0Q2;

  DilutedTrace3Factory<DilutedFactorType::Q1,
                       DilutedFactorType::Q1,
                       DilutedFactorType::Q1>
      trQ1Q1Q1;

  std::map<std::string, AbstractDilutedTraceFactory *> trace_factories;
};

class Diagram {
 public:
  Diagram(std::vector<DiagramIndex> const &corr_lookup,
          std::vector<CorrelatorRequest> const &corr_requests,
          std::string const &output_path,
          std::string const &output_filename,
          int const Lt)
      : corr_lookup_(corr_lookup),
        corr_requests_(corr_requests),
        output_path_(output_path),
        output_filename_(output_filename),
        Lt_(Lt),
        correlator_(Lt,
                    std::vector<ComplexProduct>(corr_lookup.size(), ComplexProduct{})),
        c_(omp_get_max_threads(),
           std::vector<ComplexProduct>(corr_lookup.size(), ComplexProduct{})),
        mutexes_(Lt) {}

  virtual ~Diagram() {}

  virtual char const *name() const = 0;

  std::vector<DiagramIndex> const &corr_lookup() const { return corr_lookup_; }

  std::vector<CorrelatorRequest> const &correlator_requests() const {
    return corr_requests_;
  }

  void assemble(int const t, BlockIterator const &slice_pair, DiagramParts &q) {
    int const tid = omp_get_thread_num();

    for (int i = 0; i != ssize(corr_lookup()); ++i) {
      c_[tid][i] = ComplexProduct{};
    }

    assemble_impl(c_.at(tid), slice_pair, q);

    {
      std::lock_guard<std::mutex> lock(mutexes_[t]);

      for (int i = 0; i != ssize(corr_lookup()); ++i) {
        correlator_[t][i] += c_[tid][i];
      }
    }
  }

  void write() {
    assert(output_path_ != "");
    assert(output_filename_ != "");

    WriteHDF5Correlator filehandle(
        output_path_, name(), output_filename_, make_comp_type<ComplexProduct>());

    std::vector<ComplexProduct> one_corr(Lt_);

    for (int i = 0; i != ssize(corr_lookup()); ++i) {
      for (int t = 0; t < Lt_; ++t) {
        one_corr[t] = correlator_[t][i] / static_cast<double>(Lt_);
      }
      // Write data to file.
      filehandle.write(one_corr, corr_lookup()[i]);
    }
  }

 private:
  virtual void assemble_impl(std::vector<ComplexProduct> &c,
                             BlockIterator const &slice_pair,
                             DiagramParts &q) = 0;

  std::vector<DiagramIndex> const &corr_lookup_;
  std::vector<CorrelatorRequest> const &corr_requests_;

  std::string const &output_path_;
  std::string const &output_filename_;

  int const Lt_;

  /** OpenMP-shared correlators, indices are (1) time and (2) correlator id. */
  std::vector<std::vector<ComplexProduct>> correlator_;

  /** OpenMP-shared correlators, indices are (1) thread id and (2) correlator id. */
  std::vector<std::vector<ComplexProduct>> c_;

  std::vector<std::mutex> mutexes_;
};

class GeneralDiagram : public Diagram {
 public:
  GeneralDiagram(std::vector<DiagramIndex> const &corr_lookup,
                 std::vector<CorrelatorRequest> const &corr_requests,
                 std::string const &output_path,
                 std::string const &output_filename,
                 int const Lt,
                 char const *name)
      : Diagram(corr_lookup, corr_requests, output_path, output_filename, Lt),
        name_(name) {}

  char const *name() const override { return name_; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  char const *name_;
};

/*****************************************************************************/
/*                                    C2                                     */
/*****************************************************************************/

/** Build charged 2pt correlation function
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5  \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
 *  @f}
 */
class C2c : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C2c"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

/** Build neutral 2pt correlation function
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
 *  @f}
 */
class C20 : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C20"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

/** Build neutral 2pt correlation function
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} \rangle \cdot
 *        \langle D_\mathtt{Q1}^{-1}(t'|t') \Gamma_\mathtt{Op1} \rangle
 *  @f}
 */
class C20V : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C20V"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

/*****************************************************************************/
/*                                    C3                                     */
/*****************************************************************************/

/** Build neutral 3pt correlation function
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} \rangle
 *  @f}
 */
class C3c : public Diagram {
 public:
  C3c(std::vector<DiagramIndex> const &corr_lookup,
      std::vector<CorrelatorRequest> const &corr_requests,
      std::string const &output_path,
      std::string const &output_filename,
      int const Lt);

  char const *name() const override { return "C3c"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  std::vector<std::tuple<std::array<ssize_t, 2>, std::array<ssize_t, 1>>>
      quantum_num_ids_;
};

/** Build neutral 3pt correlation function
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} \rangle
 *  @f}
 */
class C30 : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C30"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

class C30V : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C30V"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

/*****************************************************************************/
/*                                   C4                                     */
/*****************************************************************************/

/** Build charged 4pt correlation function: Direct diagram
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
 *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C4cD : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C4cD"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

/** Build neutral 4pt correlation function: Direct diagram
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
 *        \langle D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C40D : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C40D"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

/** Build charged 4pt correlation function: Vacuum diagram
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
 *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C4cV : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C4cV"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

/** Build neutral 4pt correlation function: Vacuum diagram
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
 *        \langle D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C40V : public Diagram {
 public:
  using Diagram::Diagram;

  char const *name() const override { return "C40V"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;
};

/** Build charged 4pt correlation function: Box diagram
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t'|t) \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C4cB : public Diagram {
 public:
  C4cB(std::vector<DiagramIndex> const &corr_lookup,
       std::vector<CorrelatorRequest> const &corr_requests,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4cB"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  std::vector<std::array<std::array<ssize_t, 2>, 2>> quantum_num_ids_;
};

/** Build neutral 4pt correlation function: Box diagram
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t'|t) \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C40B : public Diagram {
 public:
  C40B(std::vector<DiagramIndex> const &corr_lookup,
       std::vector<CorrelatorRequest> const &corr_requests,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C40B"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  std::vector<std::array<std::array<ssize_t, 2>, 2>> quantum_num_ids_;
};

/** Build charged 4pt correlation function: Cross diagram
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C4cC : public Diagram {
 public:
  C4cC(std::vector<DiagramIndex> const &corr_lookup,
       std::vector<CorrelatorRequest> const &corr_requests,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4cC"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  std::vector<std::array<std::array<ssize_t, 2>, 2>> quantum_num_ids_;
};

/** Build neutral 4pt correlation function: Cross diagram
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C40C : public Diagram {
 public:
  C40C(std::vector<DiagramIndex> const &corr_lookup,
       std::vector<CorrelatorRequest> const &corr_requests,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C40C"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  std::vector<std::array<std::array<ssize_t, 2>, 2>> quantum_num_ids_;
};
