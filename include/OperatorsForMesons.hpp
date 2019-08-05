#pragma once

#include "EigenVector.hpp"
#include "GaugeField.hpp"
#include "RandomVector.hpp"
#include "boost_multi_array.hpp"
#include "global_data_typedefs.hpp"
#include "typedefs.hpp"
#include "global_data.hpp"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/multi_array.hpp>

#include <Eigen/Dense>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

/** Calculates operators as they emerge in correlation functions using the
 *  stochastic estimates from the stochastic Laplacian Heaviside method
 *
 *  Basically this calculates every operator of the form
 *  - VdaggerV    : @f$ V^\dagger exp(ipx) V @f$
 *  - rVdaggerV   : @f$ (P^(b)\rho V)^\dagger exp(ipx) V @f$
 *  - rVdaggerVr  : @f$ (P^(b)\rho V)^\dagger exp(ipx) (P^(b)\rho V) @f$
 *
 *  It is straightforward to generalize the operators for Displacement, but
 *  not implemented
 *
 *  @TODO Implement derivate operators
 *  @note (MW 3.12.2017) rVdaggerVr moved to DilutedFactorFactory
 *  @note (MW 21.12.2017) rVdaggerV moved to DilutedFactorFactory
 */
class OperatorFactory {
 public:
  /** Constructor which allocates memory for all operators */
  OperatorFactory(const ssize_t Lt,
                  const ssize_t Lx,
                  const ssize_t Ly,
                  const ssize_t Lz,
                  const ssize_t nb_ev,
                  const ssize_t dilE,
                  const OperatorLookup &operator_lookuptable,
                  const std::string &handling_vdaggerv,
                  const std::string &path_vdaggerv,
                  const std::string &path_config,
                  const HypPars &hyp_parameters);
  /** Standard Destructor
   *
   *  Everything should be handled by Eigen, std::vector, and boost::multi_array
   */
  ~OperatorFactory(){};

  /****************************************************************************/
  /**************** INTERFACE FOR BUILDING ALL OPERATORS **********************/

  void create_operators(const std::string &filename,
                        const RandomVector &rnd_vec,
                        const int config,
                        const GlobalData &gd);

  void build_vdaggerv(const std::string &filename, const int config, const GlobalData &gd);
  /** Free memory of vdaggerv */
  void free_memory_vdaggerv();

  /** @todo check of vdaggerv is already build */
  inline const Eigen::MatrixXcd &return_vdaggerv(const ssize_t index,
                                                 const ssize_t t) const {
    return vdaggerv[index][t];
  }

  
 private:
  // Containers for operators which are accessible from outside
  array_Xcd_d2_eigen vdaggerv;
  /** @cond
   *  internal indices etc.
   */
  array_cd_d2 momentum;
  /** @endcond */

  /** @TODO comment private members */
  const OperatorLookup operator_lookuptable;
  const ssize_t Lt, Lx, Ly, Lz;
  const ssize_t nb_ev, dilE=0;
  bool is_vdaggerv_set = false;
  std::string handling_vdaggerv;
  std::string path_vdaggerv;
  std::string path_config;
  HypPars hyp_parameters;

  // Internal functions to build individual operators --> The interface to these
  // functions is 'create_Operators'
  // input -> filename: name and path of eigenvectors
  void read_vdaggerv(const int config);
  void read_vdaggerv_liuming(const int config);
};

void create_momenta(const ssize_t Lx,
                    const ssize_t Ly,
                    const ssize_t Lz,
                    const std::vector<VdaggerVQuantumNumbers> &vdaggerv_lookup,
                    array_cd_d2 &momentum);

void create_momenta_new(const ssize_t Lx,
                        const ssize_t Ly,
                        const ssize_t Lz,
                        const std::vector<VdaggerVQuantumNumbers> &vdaggerv_lookup,
                        array_cd_d2 &momentum);
