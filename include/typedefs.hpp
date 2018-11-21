#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/container/static_vector.hpp>
#include <boost/multi_array.hpp>

#include <algorithm>
#include <array>
#include <complex>
#include <iomanip>
#include <list>
#include <vector>

int constexpr max_rnd_ids = 10;

using RndId = int8_t;

typedef std::complex<double> Complex;

template <size_t rvecs>
using SmallVectorRndId = boost::container::static_vector<RndId, rvecs>;

template <size_t rvecs>
struct DilutedFactor {
  using Data = Eigen::MatrixXcd;

  Data data;
  std::pair<RndId, RndId> ric;
  boost::container::static_vector<RndId, rvecs> used_rnd_ids;
};

template <size_t rvecs>
struct DilutedTrace {
  using Data = Complex;

  Data data;
  boost::container::static_vector<RndId, rvecs> used_rnd_ids;
};

/** Data type for momentum */
typedef boost::multi_array<Complex, 2> array_cd_d2;
typedef boost::multi_array<Complex, 3> array_cd_d3;
typedef boost::multi_array<Complex, 4> array_cd_d4;
typedef boost::multi_array<Complex, 5> array_cd_d5;
typedef boost::multi_array<Complex, 6> array_cd_d6;
typedef boost::multi_array<Complex, 7> array_cd_d7;
typedef boost::multi_array<Complex, 8> array_cd_d8;
typedef boost::multi_array<Complex, 9> array_cd_d9;
typedef boost::multi_array<Complex, 10> array_cd_d10;

/** Special type for Quarklines */
typedef boost::multi_array<std::vector<Eigen::MatrixXcd>, 3> array_quarkline;

// Eigen typedefs
typedef std::vector<Eigen::MatrixXcd> vec_Xcd_eigen;
/** Data type for rvdaggerv and rvdaggervr */
typedef std::vector<std::vector<std::vector<Eigen::MatrixXcd>>> Xcd_d3_eigen;
/** Data type for vdaggerv */
typedef boost::multi_array<Eigen::MatrixXcd, 2> array_Xcd_d2_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 3> array_Xcd_d3_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 4> array_Xcd_d4_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 5> array_Xcd_d5_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 6> array_Xcd_d6_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 7> array_Xcd_d7_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 8> array_Xcd_d8_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 9> array_Xcd_d9_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 10> array_Xcd_d10_eigen;

/** Vector for displacing eigenvectors
 * the entries are pairs of the form (">";"x") where the first entry specifies
 * the directions forward (">") ore backward ("<") of the derivative the second
 * entry of each pair specifies the spatial direction of the derivative. One can
 * displace in "x", "y" or "z" direction.
 *
 * Example: {(>;x),(<;y),(<;x),...,(>;z)}
 */
using DisplacementDirection = std::vector<std::pair<char, char>>;

/** This is just a workaround for complex numbers to get it running for hdf5
 *
 *  @todo Change namespace into something more specific
 */
typedef struct {
  double re;
  double im;
} complex_t;

// This is the datatype to write 4pt functions and that alike directly
struct ComplexProduct {
  double rere;
  double reim;
  double imre;
  double imim;

  ComplexProduct() : rere(0.0), reim(0.0), imre(0.0), imim(0.0) {}

  ComplexProduct(const double rere,
                 const double reim,
                 const double imre,
                 const double imim)
      : rere(rere), reim(reim), imre(imre), imim(imim) {}

  ComplexProduct &operator+=(ComplexProduct const &other) {
    rere += other.rere;
    reim += other.reim;
    imre += other.imre;
    imim += other.imim;

    return *this;
  }

  ComplexProduct &operator/=(double const &right) {
    rere /= right;
    reim /= right;
    imre /= right;
    imim /= right;

    return *this;
  }

  ComplexProduct operator/(double const &right) const {
    auto res = *this;
    res /= right;
    return res;
  }
};

/** Struct to uniquely identify a sLapH operator
 *  @f$ V^\dagger exp(i(p + d/2) x) V @f$
 *
 *  In contrast to the field operator the Dirac structure is factored out
 *
 *  @todo check the expression for displacement
 *  @deprecated only displacement 0 is implement atm (26.3.17)
 */
struct VdaggerVQuantumNumbers {
  ssize_t id;
  std::array<int, 3> momentum;        /**< The -momentum as 3-vector */
  DisplacementDirection displacement; /**< The displacement as 3-vector */

  /** Constructor */
  VdaggerVQuantumNumbers(const ssize_t id,
                         std::array<int, 3> const &momentum,
                         DisplacementDirection const &displacement)
      : id(id), momentum(momentum), displacement(displacement) {}
};

/** Struct that contains all information for a sLapH operator
 *
 *  @todo confusing because rvdaggerv_lookuptable and rvdaggervr_lookuptable
 *        contain indices of vdaggerv_lookup and ric_lookup. Multiple
 *        hierarchy levels in a single struct
 *  @todo confusing because either Q1 and rvdaggerv OR Q2 and rvdaggervr are
 *        used, but in both cases half the members are spurious
 */
struct OperatorLookup {
  /** Specifies physical content of quark field operator (with Dirac structure
   *  factored out)
   */
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup;

  /** For @f$ \vec{p} = 0 @f$, @f$ V^\dagger exp(ipx) V = \mathbb{1} @f$. If
   *  applicable This contains the index of @em vdaggerv_lookup where it can
   *  be replaced by a unit matrix. If @f$ \vec{p} = 0 @f$ is not needed,
   *  @em index_of_unity is set to -1
   */
  int index_of_unity;
  bool need_gaugefield = false;

  inline ssize_t size() { return vdaggerv_lookup.size(); }
};

/** Struct that holds all information on which VdaggerV must be diluted with
 *  which random vector.
 *
 *  For rVdaggerV and rVdaggerVr the VdaggerV-operators are additionaly
 *  multiplied with random vectors. For both, VdaggerV and random index
 *  combinations there are lookuptables in OperatorLookup. This struct contains
 *  the id's of vdaggerv and ric which belong together.
 */
/** Indices needed to uniquely identify Q1 objects
 *
 *  Because vdaggerv is diagonal in dirac space, gamma may be factored out and
 *  it proves useful to calculate and reuse
 *
 *    Q1 = rvdaggerv * gamma * peram
 */
/** Indices needed to uniquely identify Q2 objects
 *
 *  Because vdaggerv is diagonal in dirac space, gamma may be factored out.
 *
 *  If the @f$ \gamma_5 @f$-trick is used it proves useful to calculate and
 *  reuse
 *
 *    Q2 = @f$ \gamma_5 @f$ peram1@f$ ^\dagger \gamma_5 @f$ * vdaggerv *
 *          gamma * peram2
 */
struct DilutedFactorIndex {
  /** Identifies physical content of @f$ V^dagger V @f$ */
  ssize_t id_vdaggerv;
  /** Flag that indicates whether VdaggerV must be daggered (prior to
   *  multiplication with random vectors) to get the correct quantum numbers
   */
  bool need_vdaggerv_daggering;

  /** List of necessarry gamma combinations */
  std::vector<int> gamma;

  /** The entries of the pair correspond to the first and second random index.
   *  List of all possible combinations of random vector indices for quarks
   *  specified by @em id_q1 and @em id_q2
   */
  std::vector<std::pair<ssize_t, ssize_t>> rnd_vec_ids;
};

inline bool operator==(DilutedFactorIndex const &first, DilutedFactorIndex const second) {
  if ((first.id_vdaggerv == second.id_vdaggerv) &&
      (first.need_vdaggerv_daggering == second.need_vdaggerv_daggering) &&
      (first.gamma == second.gamma) && (first.rnd_vec_ids == second.rnd_vec_ids))
    return true;
  else
    return false;
}

/** Maps index from CorrelatorLookup to QuarklineQ1Indicies or
 *  QuarklineQ2Indicies, depending on the quarkline needed in Correlator
 *
 *  @deprecated No longer necessary in memory_optimised branch as the quarklines
 *              are built on the fly
 */
struct DilutedFactorIndicesCollection {
  std::vector<DilutedFactorIndex> Q0;
  std::vector<DilutedFactorIndex> Q1;
  std::vector<DilutedFactorIndex> Q2V;
  std::vector<DilutedFactorIndex> Q2L;
};

// Q0 formerly called rVdaggerVr
enum class DilutedFactorType { Q0, Q1, Q2, Q2L, Q2V };

/** All information needed to build and write the correlator given the
 *  quarklines were calculated beforehand
 *
 *  - id
 *  - Indices for the quarklines and gamma structure
 *  - Paths and information for IO
 */
struct DiagramIndex {
  ssize_t id;
  std::string hdf5_dataset_name;
  std::vector<ssize_t> lookup;
  std::vector<int> gamma;
  /** Just a small constructor to ensure easy filling of its vector form */
  DiagramIndex(const ssize_t id,
               const std::string &hdf5_dataset_name,
               const std::vector<ssize_t> &lookup,
               const std::vector<int> &gamma = std::vector<int>{})
      : id(id), hdf5_dataset_name(hdf5_dataset_name), lookup(lookup), gamma(gamma){};
};

inline bool operator==(DiagramIndex const &first, DiagramIndex const &second) {
  if ((first.hdf5_dataset_name == second.hdf5_dataset_name) &&
      (first.lookup == second.lookup))
    return true;
  else
    return false;
}

template <DilutedFactorType qlt>
struct DilutedFactorTypeTraits {};

template <>
struct DilutedFactorTypeTraits<DilutedFactorType::Q0> {
  typedef std::vector<DilutedFactorIndex> type;
  static ssize_t constexpr num_times = 1;
};

template <>
struct DilutedFactorTypeTraits<DilutedFactorType::Q1> {
  typedef std::vector<DilutedFactorIndex> type;
  static ssize_t constexpr num_times = 2;
};

template <>
struct DilutedFactorTypeTraits<DilutedFactorType::Q2> {
  typedef std::vector<DilutedFactorIndex> type;
  static ssize_t constexpr num_times = 3;
};

template <>
struct DilutedFactorTypeTraits<DilutedFactorType::Q2L> {
  typedef std::vector<DilutedFactorIndex> type;
  static ssize_t constexpr num_times = 3;
};

template <>
struct DilutedFactorTypeTraits<DilutedFactorType::Q2V> {
  typedef std::vector<DilutedFactorIndex> type;
  static ssize_t constexpr num_times = 3;
};

#define MU_DEBUG(x) std::cout << std::setw(30) << #x << ": " << (x) << std::endl;

/** Special type for Correlators */
typedef boost::multi_array<std::vector<Complex>, 3> array_corr;

template <typename T>
ssize_t ssize(std::vector<T> const &vec) {
  return static_cast<ssize_t>(vec.size());
}
