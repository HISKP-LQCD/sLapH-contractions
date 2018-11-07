#include "DilutedFactorY.h"

#include <utility>
#include <omp.h>

#include <iostream>

#include "local_timer.h"

namespace {
std::complex<double> const I(0.0, 1.0);
}

template <DilutedFactorType qlt>
DilutedFactorFactory<qlt>::DilutedFactorFactory(
    RandomVector const &random_vector,
    Perambulator const &perambulator,
    OperatorFactory const &_meson_operator,
    ssize_t const dilT,
    ssize_t const dilE,
    ssize_t const nev,
    typename DilutedFactorTypeTraits<qlt>::type const &_quarkline_indices)
    : peram(perambulator),
      rnd_vec(random_vector),
      meson_operator(_meson_operator),
      dilT(dilT),
      dilE(dilE),
      nev(nev),
      quarkline_indices(_quarkline_indices) {}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// rvdaggervr is calculated by multiplying vdaggerv with the same quantum
// numbers with random vectors from right and left.
template <>
void DilutedFactorFactory<DilutedFactorType::Q0>::build(Key const &time_key) {
  int const eigenspace_dirac_size = dilD * dilE;

  auto const t1 = time_key[0];

  LT_FINE_DECLARE;

  for (int operator_key = 0; operator_key < ssize(quarkline_indices); ++operator_key) {
    auto const &op = quarkline_indices[operator_key];
    const ssize_t gamma_id = op.gamma[0];

    LT_FINE_START;
    // Using an explicit copy for vdv here makes sense because then performance
    // of the daggered and undaggered dilution is the same at the cost
    // of a single explicit transposition
    Eigen::MatrixXcd vdv;
    if( op.need_vdaggerv_daggering ){
      vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1).adjoint();
    }else{
      vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1);
    }
    LT_FINE_STOP;
    LT_FINE_PRINT( (op.need_vdaggerv_daggering ? "[DilutedFactorFactory::Q0] return_VdaggerV^dagger" :
                                                 "[DilutedFactorFactory::Q0] return_VdaggerV" ) );

    ssize_t rnd_counter = 0;
    int check = -1;
    Eigen::MatrixXcd M;  // Intermediate memory

    /*! Dilution of columns */
    for (const auto &rnd_id : op.rnd_vec_ids) {
      if (check != rnd_id.second) {  // this avoids recomputation
        LT_FINE_START;
        /*! Should be 4*nev rows, but there is always just one entry not zero */
        M = Eigen::MatrixXcd::Zero(nev, 4 * dilE);

        for (ssize_t vec_i = 0; vec_i < nev; vec_i++) {
          for (ssize_t block = 0; block < 4; block++) {
            ssize_t blk = block + (vec_i + nev * t1) * 4;
              M.col(vec_i % dilE + dilE * block).noalias() +=
                rnd_vec(rnd_id.second, blk) * vdv.col(vec_i);
          }
        }
        LT_FINE_STOP;
        LT_FINE_PRINT( (op.need_vdaggerv_daggering ? "[DilutedFactorFactory::Q0] VdaggerV^dagger*r Thread " :
                                                      "[DilutedFactorFactory::Q0] VdaggerV*r Thread ") );
      }

      LT_FINE_START;
      /*! Dilution of rows and creating a sparse matrix from smaller blocks */
      Eigen::MatrixXcd matrix =
          Eigen::MatrixXcd::Zero(eigenspace_dirac_size, eigenspace_dirac_size);
      for (ssize_t block = 0; block < 4; block++) {
        const Complex value = gamma_vec[gamma_id].value[block];
        const ssize_t gamma_index = gamma_vec[gamma_id].row[block];
        for (ssize_t vec_i = 0; vec_i < nev; vec_i++) {
          ssize_t blk = gamma_index + (vec_i + nev * t1) * dilD;
          matrix.block(vec_i % dilE + dilE * gamma_index, block * dilE, 1, dilE).noalias() +=
              value * M.block(vec_i, block * dilE, 1, dilE) *
              std::conj(rnd_vec(rnd_id.first, blk));
        }
      }
      LT_FINE_STOP;
      LT_FINE_PRINT("[DilutedFactorFactory::Q0] r*VdaggerV*r");

      check = rnd_id.second;
      rnd_counter++;
      
      Ql[time_key][{operator_key}].push_back(
          {matrix, std::make_pair(rnd_id.first, rnd_id.second), {}});
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template <>
void DilutedFactorFactory<DilutedFactorType::Q1>::build(Key const &time_key) {
  int const eigenspace_dirac_size = dilD * dilE;

  int const t1 = time_key[0];
  int const b2 = time_key[1];

  LT_FINE_DECLARE;

  for (int operator_key = 0; operator_key < ssize(quarkline_indices); ++operator_key) {
    auto const &op = quarkline_indices[operator_key];
    
    LT_FINE_START;
    // Using an explicit copy for vdv here makes sense because then performance
    // of the daggered and undaggered dilution is guaranteed to be the same
    // at the cost of a single transposition
    Eigen::MatrixXcd vdv;
    if( op.need_vdaggerv_daggering ){
      vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1).adjoint();
    }else{
      vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1);
    }
    LT_FINE_STOP;
    LT_FINE_PRINT( (op.need_vdaggerv_daggering ? 
                    "[DilutedFactorFactory::Q1] return_vdaggerv^dagger Thread " :
                    "[DilutedFactorFactory::Q1] return_vdaggerv Thread ") );

    for (auto const &rnd_id : op.rnd_vec_ids) {
      auto const gamma_id = op.gamma[0];

      LT_FINE_START;
      Eigen::MatrixXcd rvdaggerv = Eigen::MatrixXcd::Zero(eigenspace_dirac_size, nev);
      for (ssize_t vec_i = 0; vec_i < nev; ++vec_i) {
        for (ssize_t block = 0; block < dilD; block++) {
          ssize_t blk = block + vec_i * dilD + dilD * nev * t1;
            rvdaggerv.row(vec_i % dilE + dilE * block).noalias() +=
                std::conj(rnd_vec(rnd_id.first, blk)) * vdv.row(vec_i);
        }
      }
      LT_FINE_STOP;
      LT_FINE_PRINT( (op.need_vdaggerv_daggering ?
                     "[DilutedFactorFactory::Q1] r*VdaggerV^dagger Thread " :
                     "[DilutedFactorFactory::Q1] r*VdaggerV Thread ") );
      
      LT_FINE_START;
      Eigen::MatrixXcd matrix =
          Eigen::MatrixXcd::Zero(eigenspace_dirac_size, eigenspace_dirac_size);
      for (int row = 0; row < dilD; row++) {
        for (int col = 0; col < dilD; col++) {
          matrix.block(row * dilE, col * dilE, dilE, dilE).noalias() =
              gamma_vec[gamma_id].value[row] * rvdaggerv.block(row * dilE, 0, dilE, nev) *
              peram[rnd_id.second].block((t1 * dilD + gamma_vec[gamma_id].row[row]) * nev,
                                         (b2 * dilD + col) * dilE,
                                         nev,
                                         dilE);
        }
      }
      LT_FINE_STOP;
      LT_FINE_PRINT("[DilutedFactorFactory::Q1] r*VdaggerV*Peram*r");

      Ql[time_key][{operator_key}].push_back(
          {matrix, std::make_pair(rnd_id.first, rnd_id.second), {}});
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <>
void DilutedFactorFactory<DilutedFactorType::Q2>::build(Key const &time_key) {
  int const eigenspace_dirac_size = dilD * dilE;

  auto const b1 = time_key[0];
  auto const t1 = time_key[1];
  auto const b2 = time_key[2];

  LT_FINE_DECLARE;
  
  for (int operator_key = 0; operator_key < ssize(quarkline_indices); ++operator_key) {
    auto const &op = quarkline_indices[operator_key];
    ssize_t rnd_counter = 0;
    int check = -1;
    
    LT_FINE_START;
    // Using an explicit copy for vdv here makes sense because then performance
    // of the daggered and undaggered dilution is the same at the cost of a 
    // single explicit transposition
    Eigen::MatrixXcd vdv;
    if( op.need_vdaggerv_daggering ){
      vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1).adjoint();
    }else{
      vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1);
    }
    LT_FINE_STOP;
    LT_FINE_PRINT( (op.need_vdaggerv_daggering ? 
                    "[DilutedFactorFactory::Q2] return_vdaggerv^dagger Thread " :
                    "[DilutedFactorFactory::Q2] return_vdaggerv Thread ") );
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(dilD * dilE, 4 * nev);

    for (const auto &rnd_id : op.rnd_vec_ids) {
      if (check != rnd_id.first) {  // this avoids recomputation
        LT_FINE_START;
        // this is reasonably fast (~ 15 Gflop/s on a single core of lnode14)
        // but should never be computed twice because it is quite expensive
        for (int row = 0; row < dilD; row++) {
          for (int col = 0; col < 4; col++) {
              M.block(col * dilE, row * nev, dilE, nev).noalias() =
                  peram[rnd_id.first]
                      .block((t1 * 4 + row) * nev, (b1 * dilD + col) * dilE, nev, dilE)
                      .adjoint() *
                      vdv;
            // gamma_5 trick
            if (((row + col) == 3) || (abs(row - col) > 1))
              M.block(col * dilE, row * nev, dilE, nev) *= -1.;
          }
        }
        LT_FINE_STOP;
        LT_FINE_PRINT( (op.need_vdaggerv_daggering ?
                       "[DilutedFactorFactory::Q2] Peram_dagger*VdaggerV^dagger Thread " :
                       "[DilutedFactorFactory::Q2] Peram_dagger*VdaggerV Thread ") );
      }
      
      LT_FINE_START;
      Eigen::MatrixXcd matrix =
          Eigen::MatrixXcd::Zero(eigenspace_dirac_size, eigenspace_dirac_size);

      const ssize_t gamma_id = op.gamma[0];
      
      // this is reasonably fast (~15 Gflop/s) but should not be redone
      for (ssize_t block_dil = 0; block_dil < dilD; block_dil++) {
        const Complex value = gamma_vec[gamma_id].value[block_dil];
        const ssize_t gamma_index = gamma_vec[gamma_id].row[block_dil];
        for (int row = 0; row < dilD; row++) {
          for (int col = 0; col < dilD; col++) {
            matrix.block(row * dilE, col * dilE, dilE, dilE).noalias() +=
                value *
                M.block(row * dilE, block_dil * nev, dilE, nev) *
                peram[rnd_id.second].block(
                    (t1 * 4 + gamma_index) * nev, (b2 * dilD + col) * dilE, nev, dilE);
          }
        }
      }
      LT_FINE_STOP;
      LT_FINE_PRINT("[DilutedFactorFactory::Q2] r*VdaggerV*Peram*r");

      check = rnd_id.first;
      rnd_counter++;

      auto const rnd_pair = std::make_pair(rnd_id.first, rnd_id.second);

      Ql[time_key][{operator_key}].push_back({matrix, rnd_pair, {}});
    }
  }
}

template class DilutedFactorFactory<DilutedFactorType::Q0>;
template class DilutedFactorFactory<DilutedFactorType::Q1>;
template class DilutedFactorFactory<DilutedFactorType::Q2>;
