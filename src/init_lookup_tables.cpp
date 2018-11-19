/** @file
 *
 *  Functions translating lists build from the infile into lookup_tables.
 *
 *  The lookup tables only contain unique quantum number combinations and
 *  from the lists and latter are replaced by indexlists referring to those
 *  lookup tables.
 *
 *  This automatically avoids recalculation of operators or even complete
 *  correlators
 *
 *  @todo Why are most functions static and not simply in an unnamed namespace?
 */

#include "init_lookup_tables.hpp"

#include <iostream>

namespace {

using Vector = QuantumNumbers::VectorData;

/**
 *  Check whether a momentum is within the set of desired momenta specified
 *  in the correlator_list
 *
 *  @param[in]   p_tot Total 3-momentum at source or sink
 *  @param[out]  P     The specfied momenta of the reference frame
 *
 *  @returns     true if @f$|\mathtt{p_tot}|^2 \in \mathtt{P}@f$, false
 *otherwise
 */
bool desired_total_momentum(Vector const &p_tot, std::vector<Vector> const &P) {
  /** If no total momentum is specified, there is no selection.
   *  @todo It is better to force the user to specify P. Catch that in
   *        global_data_input_handling
   */
  if (P.empty()) {
    return true;
  }

  if (std::find(P.begin(), P.end(), p_tot) == P.end()) {
    return false;
  } else {
    return true;
  }
}

/**
 *  For multi-meson operators check whether the sum of momenta at the same
 *  time slice is above a given cutoff
 *
 *  @param[in]  p1  Momentum of the first meson
 *  @param[in]  p2  Momentum of the second meson
 *
 *  @returns        true if @f$|\mathtt{p1}|^2 + |\mathtt{p2}|^2 <
 *\mathtt{cutoff} @f$,
 *                  false otherwise
 *
 *  @warning Cutoff is hardcoded in this function!
 *  @warning No cutoff for P > 4 is implemented!
 */
static bool momenta_below_cutoff(Vector const &p1,
                                 Vector const &p2,
                                 std::map<int, int> const &momentum_cutoff) {
  Vector const p_tot = p1 + p2;

  if (p_tot.squaredNorm() > 4) {
    std::cout << "In momenta_below_cutoff(): WARNING! No cutoff for P > 4"
              << " implemented" << std::endl;
  }

  return p1.squaredNorm() + p2.squaredNorm() <= momentum_cutoff.at(p_tot.squaredNorm());
}

}  // end of unnamed namespace

/**
 * Build an array with all the quantum numbers needed for a particular
 * correlation function respecting physical conservation laws.
 *
 * @param[in] correlator A single correlator specified in the infile and
 * processed into the Correlators struct
 * @param[in] operator_list List of all operators specified in the infile and
 * processed into Operators struct
 * @param[out] quantum_numbers A list of all physical quantum numbers as
 * specified in the QuantumNumbers struct that are possible for @em correlator
 * @param[in] momentum_cutoff Cutoffs for sum of momenta squared, see
 * GlobalData::momentum_cutoff.
 *
 * \p correlator contains multiple operator_numbers. From combinatorics a
 * large number of combinations arise. In general only a subset of them are
 * physically allowed or necessary to calculate. In this function momentum
 * conservation is enforced and multiple cutoffs introduced.
 */
void build_quantum_numbers_from_correlator_list(
    Correlators_2 const &correlator,
    Operator_list const &operator_list,
    std::vector<std::vector<QuantumNumbers>> &quantum_numbers,
    std::map<int, int> const &momentum_cutoff) {
  std::vector<Operators> qn_op;
  for (auto const &op_number : correlator.operator_numbers) {
    if (op_number >= ssize(operator_list)) {
      std::ostringstream oss;
      oss << "Operator with ID " << op_number
          << " which is used in [correlator_lists] is not defined. Please adjust your "
             "parameter file.";
      throw std::runtime_error(oss.str());
    }
    qn_op.emplace_back(operator_list[op_number]);
  }

  std::cout << "Constructing momentum combinations for " << correlator.type << std::endl;

  // Restriction to what shall actually be computed is done in if statements
  // for each diagram because it depends on the number of quarks.
  //
  // @todo Think about a way to avoid these if conditions.
  if (correlator.type == "C1") {
    for (auto const &op0 : qn_op[0])
      quantum_numbers.emplace_back(std::vector<QuantumNumbers>({op0}));
  } else if (correlator.type == "C2c" || correlator.type == "C20" ||
             correlator.type == "C20V" || correlator.type == "Check") {
    // Build all combinations of operators and impose momentum conservation
    // and cutoffs
    for (auto const &op0 : qn_op[0]) {
      Vector p_so = op0.momentum;

      if (desired_total_momentum(p_so, correlator.tot_mom)) {
        for (const auto &op1 : qn_op[1]) {
          Vector p_si = op1.momentum;

          // momentum at source and sink must always be the same for 2pt fcts.
          if (p_so == -p_si) {
            std::vector<QuantumNumbers> single_vec_qn = {op0, op1};
            quantum_numbers.emplace_back(single_vec_qn);
          }
        }
      }
    }
  }

  else if (correlator.type == "C3c" || correlator.type == "C30") {
    std::map<int, int> counter; /** initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op2 : qn_op[2]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op2.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2, momentum_cutoff)) {
          for (const auto &op1 : qn_op[1]) {
            Vector p_si = op1.momentum;

            if (desired_total_momentum(p_si, correlator.tot_mom)) {
              if (p_so == -p_si) {
                const int p_tot = p_si.squaredNorm();
                counter[p_tot]++;

                quantum_numbers.emplace_back(std::vector<QuantumNumbers>{op0, op1, op2});
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  } else if (correlator.type == "C30V") {
    std::map<int, int> counter; /** initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op1 : qn_op[1]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op1.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2, momentum_cutoff)) {
          for (const auto &op2 : qn_op[2]) {
            Vector p_si = op2.momentum;

            if (desired_total_momentum(p_si, correlator.tot_mom)) {
              if (p_so == -p_si) {
                const int p_tot = p_si.squaredNorm();
                counter[p_tot]++;

                quantum_numbers.emplace_back(std::vector<QuantumNumbers>{op0, op1, op2});
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  }

  else if (correlator.type == "C4cD" || correlator.type == "C40D" ||
           correlator.type == "C4cC" || correlator.type == "C40C") {
    std::map<int, int> counter; /** initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op2 : qn_op[2]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op2.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2, momentum_cutoff)) {
          for (const auto &op1 : qn_op[1]) {
            for (const auto &op3 : qn_op[3]) {
              Vector p_si_1 = op1.momentum;
              Vector p_si_2 = op3.momentum;
              Vector p_si = p_si_1 + p_si_2;

              if (desired_total_momentum(p_si, correlator.tot_mom) &&
                  momenta_below_cutoff(p_si_1, p_si_2, momentum_cutoff)) {
                if (p_so == -p_si) {
                  const int p_tot = p_si.squaredNorm();
                  counter[p_tot]++;

                  // create combinations
                  quantum_numbers.emplace_back(
                      std::vector<QuantumNumbers>{op0, op1, op2, op3});
                }
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  }

  else if (correlator.type == "C4cB" || correlator.type == "C40B") {
    std::map<int, int> counter; /** initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op3 : qn_op[3]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op3.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2, momentum_cutoff)) {
          for (const auto &op1 : qn_op[1]) {
            for (const auto &op2 : qn_op[2]) {
              Vector p_si_1 = op1.momentum;
              Vector p_si_2 = op2.momentum;
              Vector p_si = p_si_1 + p_si_2;

              if (desired_total_momentum(p_si, correlator.tot_mom) &&
                  momenta_below_cutoff(p_si_1, p_si_2, momentum_cutoff)) {
                if (p_so == -p_si) {
                  const int p_tot = p_si.squaredNorm();
                  counter[p_tot]++;

                  // create combinations
                  quantum_numbers.emplace_back(
                      std::vector<QuantumNumbers>{op0, op1, op2, op3});
                }
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  }

  else if (correlator.type == "C40V" || correlator.type == "C4cV") {
    std::map<int, int> counter; /** initialized with zero */
    for (const auto &op0 : qn_op[0]) {
      for (const auto &op1 : qn_op[1]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op1.momentum;
        Vector p_so = p_so_1 + p_so_2;
        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2, momentum_cutoff)) {
          for (const auto &op2 : qn_op[2]) {
            for (const auto &op3 : qn_op[3]) {  // all combinations of operators
              Vector p_si_1 = op2.momentum;
              Vector p_si_2 = op3.momentum;
              Vector p_si = p_si_1 + p_si_2;
              if (p_so == -p_si) {
                const int p_tot = p_si.squaredNorm();
                counter[p_tot]++;

                // create combinations
                quantum_numbers.emplace_back(
                    std::vector<QuantumNumbers>{op0, op1, op2, op3});
              }
            }
          }
        }
      }
    }
    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  }
}

// TODO: Not the right place, not sure where to place it otherwise
/** Makes a string object of a displacement vector */
std::string vector_to_string(const std::vector<std::pair<char, char>> &in) {
  std::string out;
  if (in.empty())
    out = "000";
  for (auto const &dis : in) {
    out.push_back(dis.first);
    out.push_back(dis.second);
  }
  return out;
}

/** Create the names for output files and hdf5 datasets.
 *
 *  @param[in]  corr_type {C1,C2c,C20,C20V,C3c,C30,C4cD,C4cV,C4cC,C4cB,C40D,
 *                         C40V,C40C,C40B} :
 *  @param[in]  cnfg :            Number of first gauge configuration
 *  @param[in]  outpath           Output path from the infile.
 *  @param[in]  quark_types       Flavor of the quarks
 *  @param[in]  quantum_numbers   Physical quantum numbers
 *  @param[out] hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *
 *  The output path is constructed by appending a "/" to @em outpath.
 *  The output filename is built from @em corr_type and @em cnfg.
 *  The dataset name is built from @em corr_type, a letter for each
 *  @em quark_type, and the quantum numbers.
 *
 *  @todo Why don't we just build the complete path here already?
 */
static void build_correlator_names(
    const std::string &corr_type,
    int cnfg,
    const std::string &outpath,
    const std::vector<std::string> &quark_types,
    const std::vector<std::vector<QuantumNumbers>> &quantum_numbers,
    std::vector<std::string> &hdf5_dataset_name) {
  for (const auto &qn_row : quantum_numbers) {
    std::string filename = corr_type + "_";
    for (const auto &qt : quark_types)  // adding quark content
      filename += qt;
    for (const auto &qn : qn_row) {  // adding quantum numbers
      filename += std::string("_p") + to_string(qn.momentum);
      filename += std::string(".d") + to_string(qn.displacement);
      filename += std::string(".g") + to_string(qn.gamma);
    }
    hdf5_dataset_name.emplace_back(filename);
  }
}

static std::string const build_hdf5_dataset_name(
    const std::string &corr_type,
    int cnfg,
    const std::string &outpath,
    const std::vector<std::string> &quark_types,
    const std::vector<QuantumNumbers> &qn) {
  std::string filename = corr_type + "_";
  for (const auto &qt : quark_types)  // adding quark content
    filename += qt;
  for (const auto &op : qn) {  // adding quantum numbers
    filename += std::string("_p") + to_string(op.momentum);
    filename += std::string(".d") + to_string(op.displacement);
    filename += std::string(".g") + to_string(op.gamma);
  }
  return filename;
}
/******************************************************************************/
/** Translate list of QuantumNumbers into lookuptable for VdaggerV
 *
 *  @param[in]  quantum_numbers List of all quantum numbers operators are needed
 *                              for
 *  @param[out] vdaggerv_lookup Unique list of all VdaggerV operators needed.
 *                              Corresponds to @em quantum_numbers, but in
 *                              contrast does not contain Dirac structure.
 *                              Part of GlobalData::operator_lookup
 *  @param[out] vdv_indices     Indexlist referring to @em vdaggerv_lookup
 *                              to replace @em quantum_numbers
 *                              The first index is the id of VdaggerV, the
 *                              second tells us if VdaggerV must be daggered to
 *                              get the desired quantum numbers.
 */
void build_VdaggerV_lookup(
    const std::vector<std::vector<QuantumNumbers>> &quantum_numbers,
    std::vector<VdaggerVQuantumNumbers> &vdaggerv_lookup,
    std::vector<std::vector<std::pair<ssize_t, bool>>> &vdv_indices) {
  for (const auto &qn_vec : quantum_numbers) {
    std::vector<std::pair<ssize_t, bool>> vdv_indices_row;
    for (const auto &qn : qn_vec) {
      bool dagger;
      // checking if the combination of quantum numbers already exists in
      // vdaggerv. The position is stored in the iterator it.
      auto it = std::find_if(vdaggerv_lookup.begin(),
                             vdaggerv_lookup.end(),
                             [&qn, &dagger](VdaggerVQuantumNumbers vdv_qn) {
                               auto c1 = (vdv_qn.displacement == qn.displacement);
                               auto c2 = (Vector(vdv_qn.momentum.data()) == qn.momentum);
                               // also negative momentum is checked
                               auto c3 =
                                   (Vector(vdv_qn.momentum.data()) == (-1) * qn.momentum);
                               /** @TODO: Think about the daggering!! */
                               const Vector zero(0, 0, 0);
                               if (c1 and c2) {
                                 dagger = false;
                                 return true;
                               } else if ((c1 and c3) and (qn.displacement.empty())) {
                                 dagger = true;
                                 return true;
                               } else
                                 return false;
                             });
      // If the quantum number combination already exists only the id is needed
      // otherwise a new element is created at the end of the lookuptable.
      if (it != vdaggerv_lookup.end()) {
        vdv_indices_row.emplace_back((*it).id, dagger);
      } else {
        vdaggerv_lookup.emplace_back(
            VdaggerVQuantumNumbers(ssize(vdaggerv_lookup),
                                   {qn.momentum[0], qn.momentum[1], qn.momentum[2]},
                                   qn.displacement));
        vdv_indices_row.emplace_back(vdaggerv_lookup.back().id, false);
      }
    }
    vdv_indices.emplace_back(vdv_indices_row);
  }
}

/** @brief  Obtain index combinations of random vectors for charged correlator
 *          i.e. correlator utilizing @f$ \gamma_5 @f$-trick
 *
 *  @param[in]  quarks      Quarks as read from the infile and processed into
 *                          quark struct
 *  @param[in]  id_q1       Specifies which quark the first random index
 *                          belongs to
 *  @param[in]  id_q2       Specifies which quark the second random index
 *                          belongs to
 *  @param[in]  C1          Flag distinguishing whether the indexcombinations
 *                          are for C1 or not.
 *  @return                 pair of unique ids specifying quark flavor and
 *  			    number of random seed
 *
 *  For every quark propagator a statistical 1 in the form
 *  @f$ ( P^{(b)} \rho) \cdot (P^{(b)} \rho)^\dagger @f$
 *  is introduced.
 *
 *  As explained in GlobalData, when factorizing the correlators this ones
 *  are always split. To reconstruct the correct random index combinations,
 *  this function constructs all allowed combinations of random indices for
 *  a quarkline with two random indices. To avoid bias, two different random
 *  vectors must always have different seed and thus different indices.
 *
 *  The random indices are uniquely identifying quark and random vector. Thus
 *  There are @f$ \sum_i q_i N_\text{rnd}(q_i) @f$ random indices.
 *
 */
static std::vector<std::pair<ssize_t, ssize_t>> create_rnd_vec_id(
    const std::vector<quark> &quarks,
    const ssize_t id_q1,
    const ssize_t id_q2,
    const bool C1) {
  // set start and end points of rnd numbers
  auto rndq1_start = 0;
  for (auto i = 0; i < id_q1; i++)
    rndq1_start += quarks[i].number_of_rnd_vec;
  auto rndq2_start = 0;
  for (auto i = 0; i < id_q2; i++)
    rndq2_start += quarks[i].number_of_rnd_vec;

  auto rndq1_end = rndq1_start + quarks[id_q1].number_of_rnd_vec;
  auto rndq2_end = rndq2_start + quarks[id_q2].number_of_rnd_vec;

  // check if there are enough random vectors
  if ((quarks[id_q1].number_of_rnd_vec < 1) || (quarks[id_q2].number_of_rnd_vec < 1) ||
      (id_q1 == id_q2 && quarks[id_q1].number_of_rnd_vec < 2)) {
    std::cerr << "There are not enough random vectors for charged correlators"
              << std::endl;
    exit(-1);
  }

  // finally filling the array
  std::vector<std::pair<ssize_t, ssize_t>> rnd_vec_comb;
  if (!C1) {
    for (ssize_t i = rndq1_start; i < rndq1_end; ++i)
      for (ssize_t j = rndq2_start; j < rndq2_end; ++j)
        // To avoid bias, different random vectors must have different indices.
        if (i != j) {
          rnd_vec_comb.emplace_back(i, j);
        }
  } else {
    for (ssize_t i = rndq1_start; i < rndq1_end; ++i)
      // if C1 == True there is only one random vector and thus only same index
      // combinations are possible
      rnd_vec_comb.emplace_back(i, i);
  }

  return rnd_vec_comb;
}

static void build_Quarkline_lookup_one_qn(
    ssize_t const operator_id,
    std::vector<QuantumNumbers> const &quantum_numbers,
    std::vector<std::pair<ssize_t, bool>> const &vdv_indices,
    std::vector<std::pair<ssize_t, ssize_t>> const &rnd_vec_ids,
    std::vector<DilutedFactorIndex> &Ql_lookup,
    std::vector<ssize_t> &Ql_lookup_ids) {
  const auto qn = quantum_numbers[operator_id];

  const auto id_vdaggerv = vdv_indices[operator_id].first;
  const auto need_vdaggerv_daggering = vdv_indices[operator_id].second;

  // If Ql_lookup already contains the particular row and physical content,
  // just set the index to the existing QuarklineIndices, otherwise generate
  // it and set the index to the new one.
  DilutedFactorIndex const candidate{
      id_vdaggerv, need_vdaggerv_daggering, qn.gamma, rnd_vec_ids};
  auto it = std::find(Ql_lookup.begin(), Ql_lookup.end(), candidate);

  if (it != Ql_lookup.end()) {
    Ql_lookup_ids[operator_id] = it - Ql_lookup.begin();
  } else {
    Ql_lookup_ids[operator_id] = Ql_lookup.size();
    Ql_lookup.emplace_back(candidate);
  }
}

static ssize_t build_trQ1_lookup(std::vector<ssize_t> const ql_ids,
                                 std::vector<DiagramIndex> &trQ1_lookup) {
  DiagramIndex candidate(ssize(trQ1_lookup), "", ql_ids);
  auto it = std::find(trQ1_lookup.begin(), trQ1_lookup.end(), candidate);
  if (it == trQ1_lookup.end()) {
    trQ1_lookup.push_back(candidate);
    return trQ1_lookup.back().id;
  } else
    return (it - trQ1_lookup.begin());
}

/** @BUG If push_back moves the vector somewhere else, it-begin() might not
 *       give the correct id.
 */
static ssize_t build_corr0_lookup(std::vector<ssize_t> const ql_ids,
                                  std::vector<DiagramIndex> &trQ1Q1_lookup) {
  DiagramIndex candidate(ssize(trQ1Q1_lookup), "", ql_ids);
  auto it = std::find(trQ1Q1_lookup.begin(), trQ1Q1_lookup.end(), candidate);
  if (it == trQ1Q1_lookup.end()) {
    trQ1Q1_lookup.push_back(candidate);
    return trQ1Q1_lookup.back().id;
  } else
    return (it - trQ1Q1_lookup.begin());
}

static ssize_t build_corrC_lookup(std::vector<ssize_t> const ql_ids,
                                  std::vector<DiagramIndex> &trQ0Q2_lookup) {
  DiagramIndex candidate(ssize(trQ0Q2_lookup), "", ql_ids);
  auto it = std::find(trQ0Q2_lookup.begin(), trQ0Q2_lookup.end(), candidate);
  if (it == trQ0Q2_lookup.end()) {
    trQ0Q2_lookup.push_back(candidate);
    return trQ0Q2_lookup.back().id;
  } else
    return (it - trQ0Q2_lookup.begin());
}

class AbstractCandidateFactory {
 public:
  using Indices = std::vector<ssize_t>;

  AbstractCandidateFactory(std::vector<DiagramIndex> &tr_lookup,
                           Indices indices)
      : tr_lookup_(tr_lookup), indices_(indices) {}

  virtual ~AbstractCandidateFactory() {};

  virtual ssize_t make(std::vector<ssize_t> const &ql_ids) = 0;

 protected:
  std::vector<DiagramIndex> &tr_lookup_;
  Indices indices_;
};

class CandidateFactoryTrQ1 : public AbstractCandidateFactory {
 public:
  using AbstractCandidateFactory::AbstractCandidateFactory;

  ssize_t make(Indices const &ql_ids) override {
    std::vector<ssize_t> ids;
    Indices indices2;
    for (auto const index : indices_) {
      indices2.push_back(ql_ids[index]);
    }
    auto const id = build_trQ1_lookup(indices2, tr_lookup_);
    return id;
  }
};

class CandidateFactoryTrQ1Q1 : public AbstractCandidateFactory {
 public:
  using AbstractCandidateFactory::AbstractCandidateFactory;

  ssize_t make(Indices const &ql_ids) override {
    std::vector<ssize_t> ids;
    Indices indices2;
    for (auto const index : indices_) {
      indices2.push_back(ql_ids[index]);
    }
    auto const id = build_corr0_lookup(indices2, tr_lookup_);
    return id;
  }
};

// class CandidateFactoryPassthrough : public AbstractCandidateFactory {
//  public:
//   CandidateFactoryPassthrough() {}
// 
//   Indices make(Indices const &ql_ids) override {
//     return ql_ids;
//   }
// };

struct InnerLookup {
  std::vector<DilutedFactorIndex> *quarkline_lookup;
  size_t q1;
  size_t q2;
};

struct OuterLookup {
  using Factories = std::vector<std::shared_ptr<AbstractCandidateFactory>>;

  std::vector<DiagramIndex> *c_look;
  std::vector<InnerLookup> inner;
  Factories candidate_factories;
};

/**
 * Data structure containing quark lines and DilutedFactor indices.
 *
 * I really dislike the name “lookup” as as a data structure where you cannot
 * retrieve things is basically useless. But that does not mean that we can
 * name new stuff like this, I have even added it here twice!
 */
using BuildLookupLookupMap = std::map<std::string, OuterLookup>;

BuildLookupLookupMap make_build_lookup_lookup_map(GlobalData &gd) {
  BuildLookupLookupMap map;

  map["C20V"] =
      OuterLookup{&gd.correlator_lookuptable.C20V,
                  {InnerLookup{&gd.quarkline_lookuptable.Q1, 0, 0},
                   InnerLookup{&gd.quarkline_lookuptable.Q1, 1, 1}},
                  {std::shared_ptr<AbstractCandidateFactory>(new CandidateFactoryTrQ1(
                       gd.correlator_lookuptable.trQ1, std::vector<ssize_t>{0})),
                   std::shared_ptr<AbstractCandidateFactory>(new CandidateFactoryTrQ1(
                       gd.correlator_lookuptable.trQ1, std::vector<ssize_t>{1}))}};

  map["C30"] = OuterLookup{&gd.correlator_lookuptable.C30,
                           {InnerLookup{&gd.quarkline_lookuptable.Q1, 2, 0},
                            InnerLookup{&gd.quarkline_lookuptable.Q1, 0, 1},
                            InnerLookup{&gd.quarkline_lookuptable.Q1, 1, 2}},
                           {}};

  map["C3c"] = OuterLookup{&gd.correlator_lookuptable.C3c,
                           {InnerLookup{&gd.quarkline_lookuptable.Q2L, 2, 0},
                            InnerLookup{&gd.quarkline_lookuptable.Q1, 0, 1},
                            InnerLookup{&gd.quarkline_lookuptable.Q0, 1, 2}},
                           {}};

  map["C30V"] =
      OuterLookup{&gd.correlator_lookuptable.C30,
                  {InnerLookup{&gd.quarkline_lookuptable.Q1, 2, 0},
                   InnerLookup{&gd.quarkline_lookuptable.Q1, 0, 1},
                   InnerLookup{&gd.quarkline_lookuptable.Q1, 1, 2}},
                  {std::shared_ptr<AbstractCandidateFactory>(new CandidateFactoryTrQ1Q1(
                       gd.correlator_lookuptable.trQ1Q1, std::vector<ssize_t>{0, 1})),
                   std::shared_ptr<AbstractCandidateFactory>(new CandidateFactoryTrQ1(
                       gd.correlator_lookuptable.trQ1, std::vector<ssize_t>{2}))}};

  map["C40B"] = OuterLookup{&gd.correlator_lookuptable.C40B,
                            {InnerLookup{&gd.quarkline_lookuptable.Q1, 3, 0},
                             InnerLookup{&gd.quarkline_lookuptable.Q1, 0, 1},
                             InnerLookup{&gd.quarkline_lookuptable.Q1, 1, 2},
                             InnerLookup{&gd.quarkline_lookuptable.Q1, 2, 3}},
                            {}};

  map["C4cB"] = OuterLookup{&gd.correlator_lookuptable.C4cB,
                            {InnerLookup{&gd.quarkline_lookuptable.Q2L, 3, 0},
                             InnerLookup{&gd.quarkline_lookuptable.Q0, 0, 1},
                             InnerLookup{&gd.quarkline_lookuptable.Q2L, 1, 2},
                             InnerLookup{&gd.quarkline_lookuptable.Q0, 2, 3}},
                            {}};

  map["C40C"] = OuterLookup{&gd.correlator_lookuptable.C40C,
                            {InnerLookup{&gd.quarkline_lookuptable.Q1, 3, 0},
                             InnerLookup{&gd.quarkline_lookuptable.Q1, 0, 1},
                             InnerLookup{&gd.quarkline_lookuptable.Q1, 1, 2},
                             InnerLookup{&gd.quarkline_lookuptable.Q1, 2, 3}},
                            {}};

  map["C4cC"] = OuterLookup{&gd.correlator_lookuptable.C4cC,
                            {InnerLookup{&gd.quarkline_lookuptable.Q2V, 3, 0},
                             InnerLookup{&gd.quarkline_lookuptable.Q0, 0, 1},
                             InnerLookup{&gd.quarkline_lookuptable.Q2V, 1, 2},
                             InnerLookup{&gd.quarkline_lookuptable.Q0, 2, 3}},
                            {}};

  return map;
}

/**
 * Create lookuptable where to find the quarklines to build C20.
 *
 * @param[in] quarks Quarks as read from the infile and processed into quark
 * struct.
 * @param[in] quark_numbers List which quarks are specified in the infile.
 * @param[in] start_config Number of first gauge configuration.
 * @param[in] path_output Output path from the infile.
 * @param[in] quantum_numbers A list of all physical quantum numbers quantum
 * field operators for all correlators with Dirac structure factored out that
 * are possible for @em correlator.
 * @param[in] vdv_indices Indices identifying VdaggerV operators.
 * @param[out] Q1_lookup Lookuptable containing unique combinations of peram-,
 * vdv-, and ric-indices needed to built Q1.
 * @param[out] trQ1_lookup Lookuptable containign unique combinations of parts
 * tr(Q1Q1). Also known as trQ1Q1.
 * @param[out] c_look Lookup table for C20.
 */
static void build_C1_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> Q1_indices(std::vector<ssize_t>(quantum_numbers[0].size()));
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    /** Here the outgoing and incoming quarkline are identical. Thus ric_ids
     *  must get true
     */
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[0], true);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, Q1_indices);

    /** @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0
     */
    auto id1 = build_trQ1_lookup({Q1_indices[0]}, trQ1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C1", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{ssize(c_look), hdf5_dataset_name, {id1}};

    /** XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C20.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1_lookup       Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C20
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C20V_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> ql_ids(2);
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[0], true);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[1], true);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    auto id1 = build_trQ1_lookup({ql_ids[0]}, trQ1_lookup);
    auto id2 = build_trQ1_lookup({ql_ids[1]}, trQ1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C20V", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{ssize(c_look), hdf5_dataset_name, {id1, id2}};

    /** XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C2c.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q2V_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2V
 *  @param[out] trQ0Q2_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q0Q2).
 *                                Also known as trQ0Q2
 *  @param[out] c_look            Lookup table for C2c
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C2c_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q2V_lookup,
    std::vector<DiagramIndex> &trQ0Q2_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> ql_ids(2);
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    /** @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1
     */
    auto id1 = build_corrC_lookup({ql_ids[0], ql_ids[1]}, trQ0Q2_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C2c", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{ssize(c_look), hdf5_dataset_name, {id1}};

    /** XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C20.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1Q1_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C20
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C20_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> ql_ids(2);
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    /** @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1
     */
    auto id1 = build_corr0_lookup({ql_ids[0], ql_ids[1]}, trQ1Q1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C20", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{ssize(c_look), hdf5_dataset_name, {id1}};

    /** XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C30V.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1_lookup       Lookuptable containign unique combinations of
 *                                parts tr(Q1).
 *  @param[out] trQ1Q1_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C30V
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C30V_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1_lookup,
    std::vector<DiagramIndex> &trQ1Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> ql_ids(3);
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers) {
    quark_types.emplace_back(quarks[id].type);
  }

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[2], true);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    auto const id0 = build_corr0_lookup({ql_ids[0], ql_ids[1]}, trQ1Q1_lookup);
    auto const id1 = build_trQ1_lookup({ql_ids[2]}, trQ1_lookup);

    std::string const hdf5_dataset_name = build_hdf5_dataset_name(
        "C30V", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex const candidate{
        ssize(c_look), hdf5_dataset_name, {id0, id1}, std::vector<int>{}};

    /** XXX Better with std::set */
    auto const it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C30.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] c_look            Lookup table for C30
 *
 *  @bug I am fairly certain that the quarks of C30 are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_general_lookup(
    std::string const &name,
    OuterLookup const &ll,
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices) {
  std::vector<ssize_t> ql_ids(ll.inner.size());

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    for (auto const &lle : ll.inner) {
      auto const ric_ids = create_rnd_vec_id(
          quarks, quark_numbers[lle.q1], quark_numbers[lle.q2], name == "C1");
      build_Quarkline_lookup_one_qn(lle.q2,
                                    quantum_numbers[d],
                                    vdv_indices[d],
                                    ric_ids,
                                    *lle.quarkline_lookup,
                                    ql_ids);
    }

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        name, start_config, path_output, quark_types, quantum_numbers[d]);

    std::vector<ssize_t> ql_ids_new;
    if (ssize(ll.candidate_factories) == 0) {
      ql_ids_new = ql_ids;
    } else {
      for (auto const &candidate_factory : ll.candidate_factories) {
        auto const id = candidate_factory->make(ql_ids);
        ql_ids_new.push_back(id);
      }
    }

    DiagramIndex candidate(ssize(*ll.c_look), hdf5_dataset_name, ql_ids_new);

    /** XXX Better with std::set */
    auto it = std::find(ll.c_look->begin(), ll.c_look->end(), candidate);

    if (it == ll.c_look->end()) {
      ll.c_look->push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C4cD. Also sets
 *  trQ0Q2.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q2V_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2V
 *  @param[out] trQ0Q2_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q0Q2).
 *                                Also known as trQ0Q2
 *  @param[out] c_look            Lookup table for C4cD
 *
 *  C4cD like C4cC contains C2c. To reuse C2c, they all contain indices
 *  of trQ0Q2 which in turn contains the indices for rVdaggerVr and Q2.
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C4cD_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q2V_lookup,
    std::vector<DiagramIndex> &trQ0Q2_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> ql_ids(4);
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    /** @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1 / 2,3
     */
    auto id1 = build_corrC_lookup({ql_ids[0], ql_ids[1]}, trQ0Q2_lookup);
    auto id2 = build_corrC_lookup({ql_ids[2], ql_ids[3]}, trQ0Q2_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C4cD", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{ssize(c_look),
                           hdf5_dataset_name,
                           std::vector<ssize_t>({id1, id2}),
                           std::vector<int>({})};

    /** XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C40D.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1Q1_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C40D
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40D_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> ql_ids(4);
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    /** @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1 / 2,3
     */
    auto id1 = build_corr0_lookup({ql_ids[0], ql_ids[1]}, trQ1Q1_lookup);
    auto id2 = build_corr0_lookup({ql_ids[2], ql_ids[3]}, trQ1Q1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C40D", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{ssize(c_look),
                           hdf5_dataset_name,
                           std::vector<ssize_t>({id1, id2}),
                           std::vector<int>({})};

    /** XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C4cV. Also sets
 *  trQ0Q2.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q2V_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2V
 *  @param[out] trQ0Q2_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q0Q2).
 *                                Also known as trQ0Q2
 *  @param[out] c_look            Lookup table for C4cV
 *
 *  C4cV like C4cC contains C2c. To reuse C2c, they all contain indices
 *  of trQ0Q2 which in turn contains the indices for rVdaggerVr and Q2.
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C4cV_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q2V_lookup,
    std::vector<DiagramIndex> &trQ0Q2_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> ql_ids(4);
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    /** @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1 / 2,3
     */
    auto id1 = build_corrC_lookup({ql_ids[0], ql_ids[1]}, trQ0Q2_lookup);
    auto id2 = build_corrC_lookup({ql_ids[2], ql_ids[3]}, trQ0Q2_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C4cV", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{ssize(c_look),
                           hdf5_dataset_name,
                           std::vector<ssize_t>({id1, id2}),
                           std::vector<int>({})};

    /** XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/** Create lookuptable where to find the quarklines to build C40V.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1Q1_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C40V
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40V_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<ssize_t> ql_ids(4);
  std::vector<std::pair<ssize_t, ssize_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    /** @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1 / 2,3
     */
    auto id1 = build_corr0_lookup({ql_ids[0], ql_ids[1]}, trQ1Q1_lookup);
    auto id2 = build_corr0_lookup({ql_ids[2], ql_ids[3]}, trQ1Q1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C40V", start_config, path_output, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{ssize(c_look),
                           hdf5_dataset_name,
                           std::vector<ssize_t>({id1, id2}),
                           std::vector<int>({})};

    /** XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/**
 *  from GlobalData::correlators_list, GlobalData::operator::list and
 *  GlobalData::quarks
 *
 *  @bug In build_Q1_lookup the order of quarks given is consistently switched.
 */
void init_lookup_tables(GlobalData &gd) {
  for (const auto &correlator : gd.correlator_list) {
    // Build an array (quantum_numbers) with all the quantum numbers needed for
    // this particular correlation function.
    std::vector<std::vector<QuantumNumbers>> quantum_numbers;
    build_quantum_numbers_from_correlator_list(
        correlator, gd.operator_list, quantum_numbers, gd.momentum_cutoff);

    // Build the correlator and dataset names for hdf5 output files
    std::vector<std::string> quark_types;
    for (const auto &id : correlator.quark_numbers)
      quark_types.emplace_back(gd.quarks[id].type);
    std::vector<std::string> hdf5_dataset_name;
    build_correlator_names(correlator.type,
                           gd.start_config,
                           gd.path_output,
                           quark_types,
                           quantum_numbers,
                           hdf5_dataset_name);

    // Build the lookuptable for VdaggerV and return an array of indices
    // corresponding to @em quantum_numbers computed in step 1. In @em
    // vdv_indices the first entry is the id of vdv, the second tells us if vdv
    // must be daggered to get the correct quantum numbers.
    std::vector<std::vector<std::pair<ssize_t, bool>>> vdv_indices;
    build_VdaggerV_lookup(
        quantum_numbers, gd.operator_lookuptable.vdaggerv_lookup, vdv_indices);
    std::vector<std::pair<ssize_t, ssize_t>> rnd_index;

    auto const lookup_lookup_map = make_build_lookup_lookup_map(gd);

    if (correlator.type == "C30" || correlator.type == "C3c" ||
        correlator.type == "C40B" || correlator.type == "C4cB" ||
        correlator.type == "C40C" || correlator.type == "C4cC") {
      auto const &lookup_lookup = lookup_lookup_map.at(correlator.type);

      build_general_lookup(correlator.type,
                           lookup_lookup,
                           gd.quarks,
                           correlator.quark_numbers,
                           gd.start_config,
                           gd.path_output,
                           quantum_numbers,
                           vdv_indices);
    } else if (correlator.type == "C1") {
      build_C1_lookup(gd.quarks,
                      correlator.quark_numbers,
                      gd.start_config,
                      gd.path_output,
                      quantum_numbers,
                      vdv_indices,
                      gd.quarkline_lookuptable.Q1,
                      gd.correlator_lookuptable.trQ1,
                      gd.correlator_lookuptable.C1);
    } else if (correlator.type == "C2c" || correlator.type == "Check") {
      // Build the lookuptable for rVdaggerVr and return an array of indices
      // corresponding to the 'quantum_numbers' computed in step 1.
      //
      // Build the lookuptable for Q2 and return an array of indices
      // corresponding to the 'quantum_numbers' computed in step 1.
      //
      // Build the lookuptable for the correlation functions
      build_C2c_lookup(gd.quarks,
                       correlator.quark_numbers,
                       gd.start_config,
                       gd.path_output,
                       quantum_numbers,
                       vdv_indices,
                       gd.quarkline_lookuptable.Q0,
                       gd.quarkline_lookuptable.Q2V,
                       gd.correlator_lookuptable.trQ0Q2,
                       gd.correlator_lookuptable.C2c);
    } else if (correlator.type == "C4cD") {
      build_C4cD_lookup(gd.quarks,
                        correlator.quark_numbers,
                        gd.start_config,
                        gd.path_output,
                        quantum_numbers,
                        vdv_indices,
                        gd.quarkline_lookuptable.Q0,
                        gd.quarkline_lookuptable.Q2V,
                        gd.correlator_lookuptable.trQ0Q2,
                        gd.correlator_lookuptable.C4cD);
    } else if (correlator.type == "C4cV") {
      build_C4cV_lookup(gd.quarks,
                        correlator.quark_numbers,
                        gd.start_config,
                        gd.path_output,
                        quantum_numbers,
                        vdv_indices,
                        gd.quarkline_lookuptable.Q0,
                        gd.quarkline_lookuptable.Q2V,
                        gd.correlator_lookuptable.trQ0Q2,
                        gd.correlator_lookuptable.C4cV);
    } else if (correlator.type == "C20V") {
      build_C20V_lookup(gd.quarks,
                        correlator.quark_numbers,
                        gd.start_config,
                        gd.path_output,
                        quantum_numbers,
                        vdv_indices,
                        gd.quarkline_lookuptable.Q1,
                        gd.correlator_lookuptable.trQ1,
                        gd.correlator_lookuptable.C20V);
    } else if (correlator.type == "C20") {
      build_C20_lookup(gd.quarks,
                       correlator.quark_numbers,
                       gd.start_config,
                       gd.path_output,
                       quantum_numbers,
                       vdv_indices,
                       gd.quarkline_lookuptable.Q1,
                       gd.correlator_lookuptable.trQ1Q1,
                       gd.correlator_lookuptable.C20);
    } else if (correlator.type == "C30V") {
      build_C30V_lookup(gd.quarks,
                        correlator.quark_numbers,
                        gd.start_config,
                        gd.path_output,
                        quantum_numbers,
                        vdv_indices,
                        gd.quarkline_lookuptable.Q1,
                        gd.correlator_lookuptable.trQ1,
                        gd.correlator_lookuptable.trQ1Q1,
                        gd.correlator_lookuptable.C30V);
    } else if (correlator.type == "C40D") {
      build_C40D_lookup(gd.quarks,
                        correlator.quark_numbers,
                        gd.start_config,
                        gd.path_output,
                        quantum_numbers,
                        vdv_indices,
                        gd.quarkline_lookuptable.Q1,
                        gd.correlator_lookuptable.trQ1Q1,
                        gd.correlator_lookuptable.C40D);
    } else if (correlator.type == "C40V") {
      build_C40V_lookup(gd.quarks,
                        correlator.quark_numbers,
                        gd.start_config,
                        gd.path_output,
                        quantum_numbers,
                        vdv_indices,
                        gd.quarkline_lookuptable.Q1,
                        gd.correlator_lookuptable.trQ1Q1,
                        gd.correlator_lookuptable.C40V);
    } else {
      std::ostringstream oss;
      oss << "Correlator type " << correlator.type << " not known!" << std::endl;
      throw std::runtime_error(oss.str());
    }
  }

  /** Sets index_of_unity to the index of operator_lookuptable.vdaggerv_lookup
   *  where momentum and displacement are both zero, or to -1 if no such entry
   *  is found.
   */
  std::array<int, 3> const zero{0, 0, 0};
  bool found = false;
  for (const auto &op_vdv : gd.operator_lookuptable.vdaggerv_lookup)
    if ((op_vdv.momentum == zero) && (op_vdv.displacement.empty())) {
      gd.operator_lookuptable.index_of_unity = op_vdv.id;
      found = true;
    }
  if (!found)
    gd.operator_lookuptable.index_of_unity = -1;

  for (const auto &op_vdv : gd.operator_lookuptable.vdaggerv_lookup)
    if (!op_vdv.displacement.empty()) {
      gd.operator_lookuptable.need_gaugefield = true;
      break;
    }
}
