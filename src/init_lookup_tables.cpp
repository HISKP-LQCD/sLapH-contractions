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

#include "CartesianProduct.hpp"
#include "DiagramSpec.hpp"

#include <iostream>

using Vector = QuantumNumbers::VectorData;

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

  // We extract the vertices from the given diagram. As this is not critical in
  // performance we treat us with boundary checks here via `at` instead of
  // `operator[]`.
  auto const &vertices = diagram_spec.at(correlator.type).vertices;

  // The contractions that have to be done are all combinations of the various
  // possibilities from each vertex. We filter by momentum conservation later.
  // For the cartesian product of all operator sets we first need the sequence
  // of operators at each vertex.
  std::vector<ssize_t> lengths_source(vertices.first.size());
  std::transform(vertices.first.begin(),
                 vertices.first.end(),
                 lengths_source.begin(),
                 [&qn_op](int const i) { return ssize(qn_op.at(i)); });
  CartesianProduct product_source(lengths_source);

  // Same for the sink vertices.
  std::vector<ssize_t> lengths_sink(vertices.first.size());
  std::transform(vertices.first.begin(),
                 vertices.first.end(),
                 lengths_sink.begin(),
                 [&qn_op](int const i) { return ssize(qn_op.at(i)); });
  CartesianProduct product_sink(lengths_sink);

  // Space to store the quantum numbers assembled from the given index
  // combination. These are defined here to avoid many allocations.
  std::vector<QuantumNumbers> qn_source(vertices.first.size());
  std::vector<QuantumNumbers> qn_sink(vertices.second.size());
  std::vector<QuantumNumbers> qn_all(qn_source.size() + qn_sink.size());

  // We iterate over all source index combinations.
  for (auto const &indices_source : product_source) {
    // From the given index set we need to extract the operators. This is a bit
    // tricky as there are two levels of indexing here. The index `j` just
    // counts the source vertices of the diagram, always starting from 0.
    for (int j = 0; j != ssize(vertices.first); ++j) {
      // Then we need to get the vertex index (spanning source and sink) via
      // the data structure. This lets us retrieve the set of operators for
      // that vertex.
      auto const &ops = qn_op.at(vertices.first.at(j));

      // From that operator we take the one that is specified in the current
      // index set from the cartesian product.
      qn_source[j] = ops.at(indices_source[j]);

      // Perhaps slightly inconsistent the operators are not grouped by source
      // and sink in the assembly code much later. Therefore we need to copy
      // the operator into this second temporary data structure where it is
      // indexed with the global vertex index.
      qn_all[vertices.first.at(j)] = qn_source[j];
    }

    // Compute the total source momentum.
    Vector p_so{0, 0, 0};
    int sum_norm_sq = 0;
    for (auto const &q : qn_source) {
      p_so += q.momentum;
      sum_norm_sq += q.momentum.squaredNorm();
    }

    // If the operator was specified with a total momentum selection but the
    // current one is not in that list, discard the current combination of
    // operators.
    if (!correlator.tot_mom.empty() &&
        std::find(correlator.tot_mom.begin(), correlator.tot_mom.end(), p_so) ==
            correlator.tot_mom.end()) {
      continue;
    }

    // Also discard when we are beyond the momentum cutoff.
    if (sum_norm_sq > momentum_cutoff.at(p_so.squaredNorm())) {
      continue;
    }

    // The `C1` diagram is special as it is source-only. We do not need to work
    // with sink operators at all, then. Enforcing momentum conservation is not
    // sensible in this context.
    if (vertices.second.empty()) {
      quantum_numbers.push_back(qn_source);
    } else {
      // Analogously iterate through the sink vertex operators.
      for (auto const &indices_sink : product_sink) {
        for (int j = 0; j != ssize(vertices.second); ++j) {
          auto const &ops = qn_op.at(vertices.second.at(j));
          qn_sink[j] = ops.at(indices_sink[j]);
          qn_all[vertices.second.at(j)] = qn_sink[j];
        }

        Vector p_si{0, 0, 0};
        for (auto const &q : qn_sink) {
          p_si += q.momentum;
        }

        // In case that momentum is conserved we add the operators to the
        // lookup table. Here it is not grouped by source and sink but rather
        // in that global ordering.
        if (p_so == -p_si) {
          quantum_numbers.push_back(qn_all);
        }
      }
    }
  }
}

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

/**
 * Create the names for output files and hdf5 datasets.
 *
 * @param[in]  corr_type {C1,C2c,C20,C20V,C3c,C30,C4cD,C4cV,C4cC,C4cB,C40D,
 *                        C40V,C40C,C40B} :
 * @param[in]  cnfg :            Number of first gauge configuration
 * @param[in]  outpath           Output path from the infile.
 * @param[in]  quark_types       Flavor of the quarks
 * @param[in]  quantum_numbers   Physical quantum numbers
 * @param[out] hdf5_dataset_name Names for the datasets in one-to-one
 *                               correspondence to @em quantum_numbers
 *
 * The output path is constructed by appending a "/" to @em outpath.
 * The output filename is built from @em corr_type and @em cnfg.
 * The dataset name is built from @em corr_type, a letter for each
 * @em quark_type, and the quantum numbers.
 *
 * @todo Why don't we just build the complete path here already?
 */
static void build_correlator_names(
    std::string const &corr_type,
    int cnfg,
    std::string const &outpath,
    std::vector<std::string> const &quark_types,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
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
    std::string const &corr_type,
    int cnfg,
    std::string const &outpath,
    std::vector<std::string> const &quark_types,
    std::vector<QuantumNumbers> const &qn) {
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

/**
 * Translate list of QuantumNumbers into lookuptable for VdaggerV
 *
 * @param[in]  quantum_numbers List of all quantum numbers operators are needed
 *                             for
 * @param[out] vdaggerv_lookup Unique list of all VdaggerV operators needed.
 *                             Corresponds to @em quantum_numbers, but in
 *                             contrast does not contain Dirac structure.
 *                             Part of GlobalData::operator_lookup
 * @param[out] vdv_indices     Indexlist referring to @em vdaggerv_lookup
 *                             to replace @em quantum_numbers
 *                             The first index is the id of VdaggerV, the
 *                             second tells us if VdaggerV must be daggered to
 *                             get the desired quantum numbers.
 */
void build_VdaggerV_lookup(
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<VdaggerVQuantumNumbers> &vdaggerv_lookup,
    std::vector<std::vector<std::pair<ssize_t, bool>>> &vdv_indices) {
  for (auto const &qn_vec : quantum_numbers) {
    std::vector<std::pair<ssize_t, bool>> vdv_indices_row;
    for (auto const &qn : qn_vec) {
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
                               Vector const zero(0, 0, 0);
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
    std::vector<quark> const &quarks,
    ssize_t const id_q1,
    ssize_t const id_q2,
    bool const C1) {
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
  auto const qn = quantum_numbers[operator_id];

  auto const id_vdaggerv = vdv_indices[operator_id].first;
  auto const need_vdaggerv_daggering = vdv_indices[operator_id].second;

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

template <typename T>
ssize_t unique_push_back(std::vector<T> &vector, T const &element) {
  auto const it = std::find(vector.begin(), vector.end(), element);
  if (it == vector.cend()) {
    vector.push_back(element);
    return vector.cend() - 1 - vector.cbegin();
  } else {
    return it - vector.cbegin();
  }
}

struct CandidateFactory {
  std::string name;
  Indices indices;
};

ssize_t make_candidate(std::vector<Indices> &tr_lookup,
                       Indices const &indices,
                       Indices const &ql_ids) {
  std::vector<ssize_t> ids;
  Indices indices2;
  for (auto const index : indices) {
    indices2.push_back(ql_ids[index]);
  }
  auto const id = unique_push_back(tr_lookup, indices2);
  return id;
}

using Factories = std::vector<CandidateFactory>;

Factories make_candidate_factories(std::vector<std::vector<InnerLookup>> const &traces,
                                   DiagramIndicesCollection &correlator_lookuptable) {
  Factories f;

  for (auto const &trace : traces) {
    std::ostringstream name_ss;
    name_ss << "tr";

    std::vector<ssize_t> vertices;
    ssize_t previous_q2 = -1;

    for (auto const &ql : trace) {
      if (ql.name == "Q2L" || ql.name == "Q2V") {
        name_ss << "Q2";
      } else {
        name_ss << ql.name;
      }
      vertices.push_back(ql.q1);

      if (previous_q2 >= 0 && ql.q1 != previous_q2) {
        throw std::runtime_error(
            "Inconsistency in quark line definitions of the diagrams. This is not a "
            "user error but needs to be fixed by a developer.");
      }
      previous_q2 = ql.q2;
    }

    if (vertices[0] != trace.back().q2) {
      throw std::runtime_error(
          "This trace does not end with the vertex that it started with. This "
          "inconsistency needs to be fixed by a developer.");
    }

    auto const name = name_ss.str();

    if (name == "trQ1" || name == "trQ0Q2" || name == "trQ1Q1") {
      f.push_back(CandidateFactory{name, vertices});
    }
  }

  if (f.size() != 0 && f.size() != traces.size()) {
    throw std::runtime_error(
        "The number of AbstractCandidateFactory's does not match the number of "
        "traces. Either all or none of the traces must be converted.");
  }

  return f;
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
    DiagramIndicesCollection &correlator_lookuptable,
    std::map<std::string, std::vector<Indices>> &trace_indices_map,
    std::map<std::string, std::vector<DilutedFactorIndex>> &quarkline_lookup,
    OuterLookup const &ll,
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    std::string const &path_output,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<ssize_t, bool>>> const &vdv_indices) {
  size_t ql_size = 0;
  for (auto const &trace : ll.traces) {
    ql_size += trace.size();
  }
  std::vector<ssize_t> ql_ids(ql_size);
  auto &correlator_lookup = correlator_lookuptable[name];

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (auto const &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (ssize_t d = 0; d < ssize(quantum_numbers); ++d) {
    for (auto const &llx : ll.traces) {
      for (auto const &lle : llx) {
        auto const ric_ids = create_rnd_vec_id(
            quarks, quark_numbers[lle.q1], quark_numbers[lle.q2], lle.is_q1());
        build_Quarkline_lookup_one_qn(lle.q2,
                                      quantum_numbers[d],
                                      vdv_indices[d],
                                      ric_ids,
                                      quarkline_lookup[lle.name],
                                      ql_ids);
      }
    }

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        name, start_config, path_output, quark_types, quantum_numbers[d]);

    std::vector<ssize_t> ql_ids_new;
    auto const &candidate_factories = make_candidate_factories(ll.traces, correlator_lookuptable);
    if (ssize(candidate_factories) == 0) {
      ql_ids_new = ql_ids;
    } else {
      for (auto const &candidate_factory : candidate_factories) {
        auto const id = make_candidate(trace_indices_map[candidate_factory.name],
                                       candidate_factory.indices,
                                       ql_ids);
        ql_ids_new.push_back(id);
      }
    }

    DiagramIndex candidate(ssize(correlator_lookup), hdf5_dataset_name, ql_ids_new);

    unique_push_back(correlator_lookup, candidate);
  }
}

/**
 *  from GlobalData::correlators_list, GlobalData::operator::list and
 *  GlobalData::quarks
 *
 *  @bug In build_Q1_lookup the order of quarks given is consistently switched.
 */
void init_lookup_tables(GlobalData &gd) {
  for (auto const &correlator : gd.correlator_list) {
    // Build an array (quantum_numbers) with all the quantum numbers needed for
    // this particular correlation function.
    std::vector<std::vector<QuantumNumbers>> quantum_numbers;
    build_quantum_numbers_from_correlator_list(
        correlator, gd.operator_list, quantum_numbers, gd.momentum_cutoff);

    // Build the correlator and dataset names for hdf5 output files
    std::vector<std::string> quark_types;
    for (auto const &id : correlator.quark_numbers)
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

    auto const &spec = diagram_spec.at(correlator.type);
    build_general_lookup(correlator.type,
                         gd.correlator_lookuptable,
                         gd.trace_indices_map,
                         gd.quarkline_lookuptable,
                         spec,
                         gd.quarks,
                         correlator.quark_numbers,
                         gd.start_config,
                         gd.path_output,
                         quantum_numbers,
                         vdv_indices);
  }

  /** Sets index_of_unity to the index of operator_lookuptable.vdaggerv_lookup
   *  where momentum and displacement are both zero, or to -1 if no such entry
   *  is found.
   */
  std::array<int, 3> const zero{0, 0, 0};
  bool found = false;
  for (auto const &op_vdv : gd.operator_lookuptable.vdaggerv_lookup)
    if ((op_vdv.momentum == zero) && (op_vdv.displacement.empty())) {
      gd.operator_lookuptable.index_of_unity = op_vdv.id;
      found = true;
    }
  if (!found)
    gd.operator_lookuptable.index_of_unity = -1;

  for (auto const &op_vdv : gd.operator_lookuptable.vdaggerv_lookup)
    if (!op_vdv.displacement.empty()) {
      gd.operator_lookuptable.need_gaugefield = true;
      break;
    }
}
