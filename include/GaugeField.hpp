#pragma once

#include "EigenVector.hpp"
#include "RandomVector.hpp"
#include "io_utils.hpp"
#include "typedefs.hpp"

#include <boost/multi_array.hpp>
#include <unsupported/Eigen/MatrixFunctions>

#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

//! \typedef 2dim array for lookup tables
typedef boost::multi_array<int, 2> look;
//! \typedef 2dim array for gauge field matrices of one timeslice.
//!
typedef boost::multi_array<Eigen::Matrix3cd, 2> array_3cd_d2_eigen;

//! \class GaugeField
//! Class for handling the Gaugefield, displacement routines and Smearing
//! methods.
//!
//! This Class reads in timeslices from Gaugefields in lime format into a
//! Eigen::Matrix Structure. Eigenvectors or Matrices can be displaced in a
//! symmetrized and an unsymmetrized way. Furthermore Gaugetransformations can
//! be conducted either of the Gaugefield itself or passed Eigensystems.
//!
//! Smearings implemented so far are Hyp-Smearing, Ape-Smearing and Stout
//! smearing, acting directly on the read in gaugefield

class GaugeField {
 public:
  //! Constructor.

  //! t0 is first, tf is last timeslice to be read in
  GaugeField(const int Lt,
             const int Lx,
             const int Ly,
             const int Lz,
             const std::string config_path,
             const ssize_t t0,
             const ssize_t tf,
             const ssize_t ndir);
  ~GaugeField(){};
  //! Overloaded access operator.
  //!
  //! returns the complex 3x3 link variable of timeslice t at volumeindex v in
  //! direction dir
  //!
  inline const Eigen::Matrix3cd &operator()(const ssize_t t,
                                            const ssize_t v,
                                            const ssize_t dir) const {
    return (tslices.at(t))[v][dir];
  };
  //! read in gauge field
  //!
  //! slice_i is the initial, slice_f the final timeindex of the lattice.
  //! The read in lime gaugefield is internally splitted into slice_f -
  //! slice_i 3d timeslices and mapped to the eigenformat.
  //!
  void read_gauge_field(const ssize_t config_i,
                        const ssize_t slice_i,
                        const ssize_t slice_f);

  //! Hypercubic Blocking in 3 dimensions
  //!
  //! t is the timeslice index, alpha_1 the inner weight, alpha_2 the outer
  //! weight and iter the number of iterations the smearing is applied
  void smearing_hyp(const ssize_t t,
                    const double alpha_1,
                    const double alpha_2,
                    const ssize_t iter);
  //! Stout Smearing of Gaugelinks
  //!
  //! t is the timeslice index, rho is the staple weight and iter the number
  //! of applications of the smearing
  void smearing_stout(const ssize_t t, const double rho, const ssize_t iter);

  //! APE Smearing of one timeslice of Gaugelinks
  //!
  //! t is the timeslice index, alpha_1 is the staple weight, iter is the
  //! number of applications
  void smearing_ape(const ssize_t t, const double alpha_1, const ssize_t iter);

  //! Initialize lookup tables for 3d navigation through timeslices
  //!
  //! LX: Lattice extent in X direction, LY: Lattice extent in Y direction,
  //! LZ: Lattice extent in Z direction
  void init(const ssize_t LX, const ssize_t LY, const ssize_t LZ);

  // Test functions for navigation
  int get_up(const int pos, const int dir);
  int get_dn(const int pos, const int dir);

  //! Summ displacement vectors for phasefactor
  //!
  //!
  Eigen::Vector3f summed_displacement(const DisplacementDirection);

  //! Forward displacement acting to the right once in direction dir
  //!
  //! v is the eigensystem at timeslice t
  //! t is the timeslice index
  //! dir is the direction of the derivative
  Eigen::MatrixXcd forward_uv(const Eigen::VectorXcd &v,
                              const ssize_t t,
                              const int spatial_ind,
                              const int direction) const;

  //! Backward displacement acting to the right once in direction dir
  //!
  //! v is the eigensystem at timeslice t
  //! t is the timeslice index
  //! dir is the direction of the derivative
  Eigen::MatrixXcd backward_uv(const Eigen::VectorXcd &v,
                               const ssize_t t,
                               const int spatial_ind,
                               const int direction) const;

  //! Returns displaced vector or matrix
  //!
  //! v is the address of the Object to be displaced, t is the timeslice
  //! index, dir is one of 0,1 or 2 meaning x-,y- or z-direction respectively
  //! sym inidcates whether a symmetrized derivative should be used.
  // Eigen::MatrixXcd disp(const Eigen::MatrixXcd& v, const ssize_t t,
  //                      const ssize_t dir, bool sym);
  //! brief Shift eigenvectors about one step up or down in one direction
  //!
  //! v is the address of the Object to be shifted (one timeslice of an
  //! Eigensystem). step chooses if we want to shift up or down (+1 for up,-1
  //! for down), dir chooses the x,y or z direction 0,1 or 2, respectively.
  // Eigen::MatrixXcd shift(const Eigen::MatrixXcd& v, const ssize_t step,
  //                       const ssize_t dir);

  //! brief Returns displaced vector or matrix
  //
  //! v is the address of the Object to be displaced, t is the timeslice
  //! index, dir is one of 0,1 or 2 meaning x-,y- or z-direction respectively
  //! The derivative is taken in one operator $\bar{\Psi_{f}}
  //! \overleftrightarrow{D}\Psi_{f'}$.
  // Eigen::MatrixXcd symmetric_derivative(const Eigen::MatrixXcd& v,
  //                                      const ssize_t t, const ssize_t dir);
  //! brief subfunction for symmetric_derivative
  //!
  //! calculates product U_\mu(x) * V(x)
  // Eigen::MatrixXcd Umu_times_V(const Eigen::MatrixXcd& v,
  //                                              const ssize_t t,
  //                                              const ssize_t dir,
  //                                              const ssize_t verbose);
  //! brief subfunction for symmetric_derivative
  //!
  //! calculates product U_\mu(x) * V(x+\mu)
  // Eigen::MatrixXcd Umu_times_shiftedV(const Eigen::MatrixXcd& v,
  //                                              const ssize_t t,
  //                                              const ssize_t dir,
  //                                              const ssize_t verbose);

  //! Returns symmetric 2 times displaced vector or matrix
  //!
  //! v is the address of the Object to be displaced, t is the timeslice
  //! index, dir is one of 0,1 or 2 meaning x-,y- or z-direction respectively
  // Eigen::MatrixXcd disp_2(const Eigen::MatrixXcd& v, const ssize_t t,
  //                      const ssize_t dir);
  //! Gauge Transformation of timeslices
  //!
  //! For generating the transformation fields indices of the initial
  //! timeslice t0 and the final timeslice tf need to be passed
  void trafo(const ssize_t t0, const ssize_t tf);

  //! Returns a gaugetransformed transformed Eigenvector or LapH-Matrix
  //!
  //! If no transformation fields are generated, one is generated, otherwise
  //! first one is used
  Eigen::MatrixXcd trafo_ev(const Eigen::MatrixXcd &eig_sys);

  //! Returns plaquette for one timeslice
  //!
  //! t is the timeslice index of the gaugefield
  double plaque_ts(const ssize_t t);

 private:
  //! project one timeslice of lime Gaugefield to 3d timeslice
  //!
  //! t is the timeslice to be mapped timeslice points to the array of the
  //! lime Gaugefield.
  void map_timeslice_to_eigen(const ssize_t t, const double *timeslice);

  //! tmlQCD function to read in gaugefields to array of doubles
  //!
  //! gaugefield is a pointer to the storage for the gaugefield
  //! filename indicates the path to the configuration
  //! slice_i is the initial, slice_f the final timeslice to be read
  void read_lime_gauge_field_doubleprec_timeslices(double *gaugefield,
                                                   const char *filename,
                                                   const ssize_t slice_i,
                                                   const ssize_t slice_f);

  //! Build gauge array
  //!
  //! Constructs trange gaugefields, stored internally
  void build_gauge_array(const ssize_t trange);

  //! Returns palquette at one space-time point
  //!
  //! mu and nu are plaquette directions, vol is the volume index, t is the
  //! timeslice index
  //! Used internally by plaque_ts
  double plaque_pnt(const ssize_t mu,
                    const ssize_t nu,
                    const ssize_t vol,
                    const ssize_t t);

  //! Parameters specifying size of the Lattice
  const int Lt, Lx, Ly, Lz, V3, dim_row, V_TS, V_for_lime;
  //! Path to gauge configurations
  const std::string config_path;
  //! Vector holding boost multi_arrays for range of timeslices
  std::vector<array_3cd_d2_eigen> tslices;
  //! One 3d timeslice of Gaugelinks as 2d boost multi_array of 3x3
  //! Matrices
  array_3cd_d2_eigen omega;
  //! 2d boost_multiarray for indices in up direction
  look iup;

  //! 2d boost_multiarray for indices in down direction
  look idown;
};
