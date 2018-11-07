#include "GaugeField.hpp"

#include <boost/format.hpp>

#include <sstream>
#include <vector>

// member initializer is executed from left to right. Is used to set constant members
GaugeField::GaugeField(const int _Lt,
                       const int _Lx,
                       const int _Ly,
                       const int _Lz,
                       const std::string _config_path,
                       const ssize_t t0,
                       const ssize_t tf,
                       const ssize_t ndir)
    : Lt(_Lt),
      Lx(_Lx),
      Ly(_Ly),
      Lz(_Lz),
      V3(Lx * Ly * Lz),
      dim_row(V3 * 3),
      V_TS(dim_row * 4 * 3 * 2),
      V_for_lime(V_TS * Lt),
      config_path(_config_path),
      tslices(),
      iup(),
      idown() {
  iup.resize(boost::extents[V3][ndir]);
  idown.resize(boost::extents[V3][ndir]);
  init(Lx, Ly, Lz);
  tslices.resize(tf - t0 + 1);
  for (auto &t : tslices)
    t.resize(boost::extents[V3][ndir]);
}

// Initialize the lookup tables

void GaugeField::init(const ssize_t L1, const ssize_t L2, const ssize_t L3) {
  int *x0_h = new int[3];
  int *x1_h = new int[3];
  int *x2_h = new int[3];

  int L0_h = L1 * L2;

  for (int x0 = 0; x0 < L1; ++x0) {  // loop x0
    x0_h[2] = x0 * L0_h;
    // negative direction (index at lower boundary)
    if ((x0_h[0] = x0 - 1) < 0)
      x0_h[0] = L0_h * (L1 - 1);
    else
      x0_h[0] *= L0_h;
    // positive direction (index at upper boundary)
    if ((x0_h[1] = x0 + 1) == L1)
      x0_h[1] = 0;
    else
      x0_h[1] *= L0_h;

    for (int x1 = 0; x1 < L1; ++x1) {  // loop x1
      x1_h[2] = x1 * L2;
      // neg. dir.
      if ((x1_h[0] = x1 - 1) < 0)
        x1_h[0] = L2 * (L1 - 1);
      else
        x1_h[0] *= L2;
      // pos. dir.
      if ((x1_h[1] = x1 + 1) == L1)
        x1_h[1] = 0;
      else
        x1_h[1] *= L2;

      for (int x2 = 0; x2 < L2; ++x2) {  // loop x2
        x2_h[2] = x2;
        // neg. dir.
        if ((x2_h[0] = x2 - 1) < 0)
          x2_h[0] = L2 - 1;
        // pos. dir.
        if ((x2_h[1] = x2 + 1) == L2)
          x2_h[1] = 0;
        // overall volume index
        int i = x0_h[2] + x1_h[2] + x2_h[2];
        // std::cout << x0 << " " << x1 << " " << x2 << " " << i << std::endl;
        // upwards
        iup[i][0] = x0_h[1] + x1_h[2] + x2_h[2];
        iup[i][1] = x0_h[2] + x1_h[1] + x2_h[2];
        iup[i][2] = x0_h[2] + x1_h[2] + x2_h[1];
        // downwards
        idown[i][0] = x0_h[0] + x1_h[2] + x2_h[2];
        idown[i][1] = x0_h[2] + x1_h[0] + x2_h[2];
        idown[i][2] = x0_h[2] + x1_h[2] + x2_h[0];
      }  // end loop x2
    }    // end loop x1
  }      // end loop x0
  delete x0_h;
  delete x1_h;
  delete x2_h;
}

// get index in positive direction
int GaugeField::get_up(const int pos, const int dir) {
  return iup[pos][dir];
}
// get index in negative direction
int GaugeField::get_dn(const int pos, const int dir) {
  return idown[pos][dir];
}

// Transform one timeslice from lime array to array_3cd_d2_eigen

// mapping from gauge config to Eigen 3x3 complex matrix arrays
void GaugeField::map_timeslice_to_eigen(const ssize_t t, const double *timeslice) {
  // Number of directions
  const int NDIR = 4;
  // Number of colors
  const int NCOL = 3;

  // read in elements
  int el_input = 0;
  for (int z = 0; z < Lz; ++z) {  // spatial loops
    for (int y = 0; y < Ly; ++y) {
      for (int x = 0; x < Lx; ++x) {
        for (int mu = 1; mu < 4; ++mu) {  // direction loop
          std::complex<double> array[9];
          for (int a = 0; a < 3; ++a) {  // colour loops
            for (int b = 0; b < 3; ++b) {
              // timeslice index of real part
              int ind_r = z * V_TS / Lz + y * V_TS / (Lz * Ly) + x * V_TS / (V3) +
                          mu * V_TS / (V3 * NDIR) + a * V_TS / (V3 * NDIR * NCOL) +
                          b * V_TS / (V3 * NDIR * NCOL * NCOL) + 0;
              // timeslice index of imaginary part
              int ind_i = z * V_TS / Lz + y * V_TS / (Lz * Ly) + x * V_TS / (V3) +
                          mu * V_TS / (V3 * NDIR) + a * V_TS / (V3 * NDIR * NCOL) +
                          b * V_TS / (V3 * NDIR * NCOL * NCOL) + 1;

              std::complex<double> pair(timeslice[ind_r], timeslice[ind_i]);
              // array to be mapped to Eigen Array
              array[3 * b + a] = pair;
              ++el_input;
            }
          }
          Eigen::Map<Eigen::Matrix3cd> dummy(array);
          // spatial index
          int ind = z * Ly * Lx + y * Lx + x;
          tslices.at(t)[ind][mu - 1] = dummy;
        }
      }
    }
  }
}

// helpers for smearing and gaugetrafos

std::array<int, 2> get_dirs(int mu) {
  std::array<int, 2> dirs;
  if (mu == 0) {
    dirs.at(0) = 1;
    dirs.at(1) = 2;
  } else if (mu == 1) {
    dirs.at(0) = 0;
    dirs.at(1) = 2;
  } else if (mu == 2) {
    dirs.at(0) = 0;
    dirs.at(1) = 1;
  }

  return (dirs);
}

// storage map for dec_timeslice
static int decor(int dir, int smear_plane) {
  int ret;
  if (dir == 0) {
    if (smear_plane == 1)
      ret = 0;
    else
      ret = 3;
  } else if (dir == 1) {
    if (smear_plane == 0)
      ret = 1;
    else
      ret = 4;
  } else {
    if (smear_plane == 0)
      ret = 2;
    else
      ret = 5;
  }
  return ret;
}

// checks if su is SU(3) matrix by tr(x*x^dagger) and sum(x*x^dagger)
static void check_su3(const Eigen::Matrix3cd &su) {
  std::complex<double> test1 = (su * su.adjoint()).trace();
  std::complex<double> test2 = (su * su.adjoint()).sum();

  if (((abs(test1) - 3.) > 10e-11) || ((abs(test2) - 3.) > 10e-11)) {
    std::cout << std::setprecision(13) << " det:  " << su.determinant()
              << "\n tr: " << test1 << "\n sum: " << test2 << "\n"
              << std::endl;
  }
}

static Eigen::Matrix3cd proj_to_su3_imp(Eigen::Matrix3cd &in) {
  // avoid possible aliasing issues:
  Eigen::Matrix3cd out = ((in.adjoint() * in).sqrt()).inverse();
  in = in * out;
  std::complex<double> det = 1. / pow(in.determinant(), 1. / 3.);
  in *= det;
  return in;
}
// Build a SU(3)-matrix from a random complex matrix and a back-projection
static Eigen::Matrix3cd construct_random_su3() {
  std::vector<std::complex<double>> random(9);
  // 9 complex random numbers
  for (auto &c : random)
    c = std::complex<double>((double(rand()) / RAND_MAX) - 0.5,
                             (double(rand()) / RAND_MAX) - 0.5);
  Eigen::Matrix3cd y =
      Eigen::Matrix<std::complex<double>, 3, 3, Eigen::RowMajor>(random.data());
  y = proj_to_su3_imp(y);
  check_su3(y);
  return y;
}

// Smearing methods

// Stout-Smearing
void GaugeField::smearing_stout(const ssize_t t, const double rho, const ssize_t iter) {
  std::complex<double> im_half(0, 0.5);
  array_3cd_d2_eigen eigen_timeslice_ts(boost::extents[V3][3]);
  // Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3];
  // for ( auto i = 0; i < V3; ++i ) {
  //   eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
  // }
  for (ssize_t j = 0; j < iter; ++j) {
    for (ssize_t i = 0; i < V3; ++i) {
      for (ssize_t dir = 0; dir < 3; ++dir) {
        Eigen::Matrix3cd staple =
            Eigen::Matrix3cd::Zero();  // Holding all smearing matrices for one link
        // filling each element of smearer using 3 links
        // For each link calculate staples summing them up
        for (ssize_t not_dir = 0; not_dir < 3; ++not_dir) {
          if (dir != not_dir) {
            int mu = iup[i][not_dir];
            int nu = iup[mu][dir];
            int eta = idown[nu][not_dir];
            // Staples in positive direction
            staple +=
                (tslices.at(t))[i][not_dir] *
                ((tslices.at(t))[mu][dir] * ((tslices.at(t))[eta][not_dir].adjoint()));

            mu = idown[i][not_dir];
            nu = iup[mu][dir];
            // Staples in negative direction
            staple += ((tslices.at(t))[mu][not_dir].adjoint()) *
                      ((tslices.at(t))[mu][dir] * (tslices.at(t))[nu][dir]);
          }
        }
        Eigen::Matrix3cd omega =
            (staple * (rho / 4.)) * ((tslices.at(t))[i][dir].adjoint());
        Eigen::Matrix3cd q =
            ((omega.adjoint() - omega) * 0.5) -
            ((((omega.adjoint() - omega).trace()) * Eigen::Matrix3cd::Identity()) *
             (1. / 6.));
        eigen_timeslice_ts[i][dir] = ((q * (-1)).exp()) * (tslices.at(t))[i][dir];
      }
    }
    for (ssize_t i = 0; i < V3; ++i) {
      for (ssize_t mu = 0; mu < 3; ++mu) {
        tslices.at(t)[i][mu] = eigen_timeslice_ts[i][mu];
      }
    }
  }
}

// APE-Smearing
void GaugeField::smearing_ape(const ssize_t t, const double alpha_1, const ssize_t iter) {
  // Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3];
  // for ( auto i = 0; i < V3; ++i ) {
  //  eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
  //}
  // temporal timeslice from decorated links
  array_3cd_d2_eigen eigen_timeslice_ts(boost::extents[V3][3]);
  for (ssize_t j = 0; j < iter; ++j) {
    for (ssize_t i = 0; i < V3; ++i) {
      for (ssize_t dir = 0; dir < 3; ++dir) {
        Eigen::Matrix3cd staple =
            Eigen::Matrix3cd::Zero();  // Holding all smearing matrices for one link
        // filling each element of smearer using 3 links
        // Position indices mu nu eta
        for (ssize_t not_dir = 0; not_dir < 3; ++not_dir) {
          if (dir != not_dir) {
            int mu = iup[i][not_dir];
            int nu = iup[mu][dir];
            int eta = idown[nu][not_dir];
            // Staples in positive direction
            staple +=
                (tslices.at(t))[i][not_dir] *
                ((tslices.at(t))[mu][dir] * ((tslices.at(t))[eta][not_dir].adjoint()));

            mu = idown[i][not_dir];
            nu = iup[mu][dir];
            // Staples in negative direction
            staple += (tslices.at(t)[mu][not_dir].adjoint()) *
                      ((tslices.at(t))[mu][dir] * (tslices.at(t))[nu][dir]);
          }
        }
        eigen_timeslice_ts[i][dir] =
            (tslices.at(t)[i][dir] * (1. - alpha_1)) + (staple * alpha_1 / 4.);
      }
    }
    for (ssize_t i = 0; i < V3; ++i) {
      for (ssize_t mu = 0; mu < 3; ++mu) {
        tslices.at(t)[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
      }
    }
  }
}

// HYP-Smearing

void GaugeField::smearing_hyp(const ssize_t t,
                              const double alpha_1,
                              const double alpha_2,
                              const ssize_t iter) {
  // temporal timeslice twice the size for decorated links
  array_3cd_d2_eigen dec_timeslice(boost::extents[V3][6]);

  // temporal timeslice from decorated links
  array_3cd_d2_eigen eigen_timeslice_ts(boost::extents[V3][3]);
  // temporal integers holding directions for decorated smearing
  int mu, nu, eta;
  for (ssize_t run = 0; run < iter; ++run) {
    // calculate inner staple from original timeslice, store in dec_timeslice each link
    // can get smeared in two planes
    for (int vol = 0; vol < V3; ++vol) {
      for (int dir = 0; dir < 3; ++dir) {
        // inner staple
        Eigen::Matrix3cd inner_staple = Eigen::Matrix3cd::Zero();
        std::array<Eigen::Matrix3cd, 2> tmp_staples;
        std::array<int, 2> perpendics = get_dirs(dir);
        for (auto it_perp_dir = perpendics.begin(); it_perp_dir != perpendics.end();
             ++it_perp_dir) {
          int perp_dir = *it_perp_dir;
          // up-type smearing_indices
          mu = iup[vol][perp_dir];
          nu = iup[mu][dir];
          eta = iup[vol][dir];
          // product of up matrices
          inner_staple =
              tslices.at(t)[vol][perp_dir] *
              (tslices.at(t)[mu][dir] * (tslices.at(t)[eta][perp_dir].adjoint()));
          // down-type smearing indices
          mu = idown[vol][perp_dir];
          nu = iup[mu][dir];
          // eta is same endpoint no adjoint necessary here
          // product of down matrices
          inner_staple += (tslices.at(t)[mu][perp_dir].adjoint()) *
                          (tslices.at(t)[mu][dir] * tslices.at(t)[nu][perp_dir]);
          // Careful placement of decorated links in dec_timeslices:
          // dir=0 has placement in dec_dir = 0 (smeared in 1 plane)
          //                       dec_dir = 3 (smeared in 2 plane)
          // dir=1 has placement in dec_dir = 1 (smeared in 0 plane)
          //                       dec_dir = 4 (smeared in 2 plane)
          // dir=2 has placement in dec_dir = 2 (smeared in 0 plane)
          //                       dec_dir = 5 (smeared in 1 plane)
          Eigen::Matrix3cd stac =
              (tslices.at(t)[vol][dir] * (1 - alpha_2)) + (inner_staple * alpha_2 / 2.);
          int n_el = it_perp_dir - perpendics.begin();
          tmp_staples.at(n_el) = proj_to_su3_imp(stac);
        }

        // staple link in direction dir in non participating and negative directions
        dec_timeslice[vol][dir] = tmp_staples.at(0);
        dec_timeslice[vol][dir + 3] = tmp_staples.at(1);
      }
    }
    // calculate outer staple from dec_timeslice as modified ape-smearing

    for (int i = 0; i < V3; ++i) {
      for (int dir = 0; dir < 3; ++dir) {
        Eigen::Matrix3cd outer_staple =
            Eigen::Matrix3cd::Zero();  // Holding all smearing matrices for one link

        // filling each element of smearer using 3 links
        // debugging link
        Eigen::Matrix3cd outer_staple_test;
        for (int not_dir = 0; not_dir < 3; ++not_dir) {
          if (dir != not_dir) {
            // calculate plane in which was smeared
            int plane = ((dir + 1) ^ (not_dir + 1)) - 1;
            mu = iup[i][not_dir];
            nu = iup[mu][dir];
            eta = idown[nu][not_dir];

            // Staples in positive direction
            // replace directions by appropriate decor
            int a, b;
            a = decor(not_dir, plane);
            b = decor(dir, plane);
            outer_staple += dec_timeslice[i][a] *
                            (dec_timeslice[mu][b] * (dec_timeslice[eta][a].adjoint()));
            mu = idown[i][not_dir];
            nu = iup[mu][dir];

            // Staples in negative direction
            outer_staple += (dec_timeslice[mu][a].adjoint()) *
                            (dec_timeslice[mu][b] *
                             dec_timeslice[nu][a]);  // structure has to be a_dag, b, a
          }
        }
        eigen_timeslice_ts[i][dir] =
            (tslices.at(t)[i][dir] * (1. - alpha_1)) + (outer_staple * alpha_1 / 4.);
      }
    }
    for (int i = 0; i < V3; ++i) {
      for (int mu = 0; mu < 3; ++mu) {
        tslices.at(t)[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
      }
    }
  }
  // clean up
}

// Displacement routines, returning one Eigenvector/Eigensystem

/*
//Derivative, toogle symmetrization via sym
Eigen::MatrixXcd GaugeField::disp(const Eigen::MatrixXcd& v,
                                     const ssize_t t, const ssize_t dir, bool forward ) {

  //Information on Matrix size
  const int dim_col = v.cols();
  //Loop over all eigenvectors in
    //storing eigenvector
    Eigen::VectorXcd in(dim_row);
    Eigen::MatrixXcd out(dim_row, dim_col);
  for(int ev=0; ev < dim_col; ++ev){
    in = v.col(ev);

    //Displace eigenvector
    for (int spatial_ind = 0; spatial_ind < V3; ++spatial_ind) {
      //std::cout << "x: " << spatial_ind << std::endl;
      Eigen::Vector3cd tmp;
      Eigen::Vector3cd quark_up;
      Eigen::Vector3cd quark_down;

      //determine needed indices from lookup tables;
      int up_ind = iup[spatial_ind][dir];
      int down_ind = idown[spatial_ind][dir];

      quark_up = in.segment(3*up_ind,3);
      quark_down = in.segment(3*down_ind,3);
      if(forward) {
        tmp = 0.5 * ( ( (tslices.at(t))[spatial_ind][dir] * quark_up) -
            ( ( (tslices.at(t))[down_ind][dir].adjoint() ) * quark_down) );
      }
      // Do we need to introduce the backward derivative?
      else {
        tmp = 0.5 * ( ( (tslices.at(t))[spatial_ind][dir] * quark_up) -
            ( ( (tslices.at(t))[down_ind][dir].adjoint() ) * quark_down) );
        //Eigen::Vector3cd quark_point = in.segment(3*spatial_ind,3);
        //tmp = ( (tslices.at(t))[spatial_ind][dir] * quark_up) - quark_point ;
      }
      (out.col(ev)).segment(3*spatial_ind,3) = tmp;
    }//end spatial loop
  }//end eigenvector loop
  return out;
}

Eigen::MatrixXcd GaugeField::shift(const Eigen::MatrixXcd& v,
                                   const ssize_t step,
                                   const ssize_t dir){
  look lookuptable;
  // +1 means we want to shift up
  if(step == 1) lookuptable = iup;
  // -1 means we want to shift down
  if(step == -1) lookuptable = idown;
  //Information on Matrix size
  const int dim_col = v.cols();
  const int dim_row = v.rows();
  Eigen::MatrixXcd out(dim_row,dim_col);
  for(int ev=0; ev < dim_col; ++ev){
    for (int spatial_ind = 0; spatial_ind < V3; ++spatial_ind) {
      int lookup_ind = lookuptable[spatial_ind][dir];
      (out.col(ev)).segment(3*spatial_ind,3) = (v.col(ev)).segment(3*lookup_ind,3);
    }
  }
  return out;
}

// Symmetric Derivative
// Constructs 0.5*( V^\dag(x)U_mu(x)V(x+\mu)
//                - V^\dag(x-mu)U^\dag_mu(x-\mu)V(x) )
//
Eigen::MatrixXcd GaugeField::symmetric_derivative(const Eigen::MatrixXcd& v,
                                                  const ssize_t t,
                                                  const ssize_t dir) {
  //Information on Matrix size
  const int dim_col = v.cols();
  const int dim_row = v.rows();
  //Loop over all eigenvectors in
  //storing eigenvector
  Eigen::MatrixXcd out(dim_row, dim_col);
  // TODO: Should come out without that?
  Eigen::MatrixXcd w_dagger(dim_col, dim_row);
  //Shift goes in negative direction -> step == -1
  w_dagger = shift(Umu_times_V(v,t,dir,0),-1,dir).adjoint();
  out = 0.5*(v.adjoint()*Umu_times_shiftedV(v,t,dir,0) - w_dagger*v);
  return out;
}
*/

static Eigen::Vector3f cartesian_vector(const char dim) {
  Eigen::Vector3f cart(0, 0, 0);
  switch (dim) {
    case 'x':
      cart = Eigen::Vector3f::UnitX();
      break;
    case 'y':
      cart = Eigen::Vector3f::UnitY();
      break;
    case 'z':
      cart = Eigen::Vector3f::UnitZ();
      break;
      throw std::runtime_error("Dimension unknown");
  }
  return cart;
}

Eigen::Vector3f GaugeField::summed_displacement(
    const DisplacementDirection displacement) {
  Eigen::Vector3f k;
  k = Eigen::Vector3f::Zero();
  for (const auto &d : displacement) {
    // TODO: Think that this does not work, d.second is of type char
    d.first == '>' ? k += cartesian_vector(d.second) : k -= cartesian_vector(d.second);
  }
  return k;
}

// Calculates U_mu(x)V(x+\hat{\mu})
Eigen::MatrixXcd GaugeField::forward_uv(const Eigen::VectorXcd &v,
                                        const ssize_t t,
                                        const int spatial_ind,
                                        const int direction) const {
  int up_ind = iup[spatial_ind][direction];

  return (tslices.at(t))[spatial_ind][direction] * v.segment(3 * up_ind, 3);
}

// Calculates U_mu^dagger(x-\hat{mu})V(x-\hat{\mu})
Eigen::MatrixXcd GaugeField::backward_uv(const Eigen::VectorXcd &v,
                                         const ssize_t t,
                                         const int spatial_ind,
                                         const int direction) const {
  int down_ind = idown[spatial_ind][direction];

  return (tslices.at(t))[down_ind][direction].adjoint() * v.segment(3 * down_ind, 3);
}

/*
// Calculate U_mu(x)V(x)
Eigen::MatrixXcd GaugeField::Umu_times_V(const Eigen::MatrixXcd& v,
                                                  const ssize_t t,
                                                  const ssize_t dir,
                                                  const ssize_t verbose) {
  //Information on Matrix size
  const int dim_col = v.cols();
  const int dim_row = v.rows();
  //Loop over all eigenvectors in
  //storing eigenvector
  Eigen::VectorXcd in(dim_row);
  Eigen::MatrixXcd out(dim_row, dim_col);
  if (verbose){
    std::cout << t<< std::endl;
    std::cout << in.rows() << std::endl;
    std::cout << out.rows() << std::endl;
  }
  for(int ev=0; ev < dim_col; ++ev){
    in = v.col(ev);
    //multiply eigenvector with according gauge matrix
    for (int spatial_ind = 0; spatial_ind < V3; ++spatial_ind) {
      (out.col(ev)).segment(3*spatial_ind,3)=(tslices.at(t))[spatial_ind][dir]
                                              *in.segment(3*spatial_ind,3);
    }//end spatial loop
  }//end eigenvector loop
  return out;
}

// Calculates U_mu(x)V(x+\hat{\mu})
Eigen::MatrixXcd GaugeField::Umu_times_shiftedV(const Eigen::MatrixXcd& v,
                                                  const ssize_t t,
                                                  const ssize_t dir,
                                                  const ssize_t verbose) {
  //Information on Matrix size
  const int dim_col = v.cols();
  const int dim_row = v.rows();
  //Loop over all eigenvectors in
  //storing eigenvector
  Eigen::VectorXcd in(dim_row);
  Eigen::MatrixXcd out(dim_row, dim_col);
  for(int ev=0; ev < dim_col; ++ev){
    in = v.col(ev);
    //multiply eigenvector with according gauge matrix
    for (int spatial_ind = 0; spatial_ind < V3; ++spatial_ind) {
      Eigen::Vector3cd quark_up;
      int up_ind = iup[spatial_ind][dir];
      quark_up = in.segment(3*up_ind,3);
      (out.col(ev)).segment(3*spatial_ind,3)=(tslices.at(t))[spatial_ind][dir]
                                              *quark_up;
    }//end spatial loop
  }//end eigenvector loop
  return out;
}

///////////////////////////////////////////////////////////////////////////////
///Liumings second derivative returning one Eigenvector/Eigensystem////////////
///////////////////////////////////////////////////////////////////////////////
//Derivative, toogle symmetrization via sym
Eigen::MatrixXcd GaugeField::disp_2(const Eigen::MatrixXcd& v,
                                     const ssize_t t, const ssize_t dir) {

  //Information on Matrix size
  const int dim_col = v.cols();
  //Loop over all eigenvectors in
    //storing eigenvector
    Eigen::VectorXcd in(dim_row);
    Eigen::MatrixXcd out(dim_row, dim_col);
  for(int ev=0; ev < dim_col; ++ev){
    in = v.col(ev);

    //Displace eigenvector
    for (int spatial_ind = 0; spatial_ind < V3; ++spatial_ind) {
      //std::cout << "x: " << spatial_ind << std::endl;
      Eigen::Vector3cd tmp;
      Eigen::Vector3cd quark_up;
      Eigen::Vector3cd quark_down;
      Eigen::Vector3cd quark_double_up;
      Eigen::Vector3cd quark_double_down;

      //determine needed indices from lookup tables;
      int up_ind = iup[spatial_ind][dir];
      int down_ind = idown[spatial_ind][dir];
      int double_up_ind = iup[up_ind][dir];
      int double_down_ind = idown[down_ind][dir];

      quark_up = in.segment(3*up_ind,3);
      quark_down = in.segment(3*down_ind,3);
      quark_double_up = in.segment(3*double_up_ind,3);
      quark_double_down = in.segment(3*double_down_ind,3);
      // Holds only for symmetric derivatives
      tmp =
0.25*(tslices.at(t)[spatial_ind][dir]*tslices.at(t)[up_ind][dir]*quark_double_up
                  -(tslices.at(t)[down_ind][dir].adjoint()*tslices.at(t)[down_ind][dir]
                  +tslices.at(t)[spatial_ind][dir]*tslices.at(t)[spatial_ind][dir].adjoint())*in.segment(3*spatial_ind,3)
                  +tslices.at(t)[down_ind][dir].adjoint()*tslices.at(t)[double_down_ind][dir].adjoint()*in.segment(3*double_down_ind,3));
      (out.col(ev)).segment(3*spatial_ind,3) = tmp;
    }//end spatial loop
  }//end eigenvector loop
  return out;
}
*/

// Gaugefield transformations

// build gauge array
void GaugeField::build_gauge_array(const ssize_t trange) {
  // parameter passing still to be improved
  const ssize_t V3 = Lx * Ly * Lz;
  srand(1227);
  // resize omega
  omega.resize(boost::extents[trange][V3]);
  // fill omega
  for (ssize_t t = 0; t < trange; ++t) {
    for (ssize_t vol = 0; vol < V3; ++vol) {
      omega[t][vol] = construct_random_su3();
    }
  }
}
// Gauge-transform config for every timeslice
void GaugeField::trafo(const ssize_t t0, const ssize_t tf) {
  const ssize_t V3 = Lx * Ly * Lz;

  build_gauge_array(tf - t0);
  for (ssize_t t = t0; t < tf; ++t) {
    for (ssize_t v = 0; v < V3; ++v) {
      for (ssize_t mu = 0; mu < 3; ++mu) {
        int w = iup[v][mu];
        tslices.at(t)[v][mu] =
            omega[t - t0][v].adjoint() * (tslices.at(t)[v][mu] * omega[t - t0][w]);
      }
    }
  }
}

// TODO: work on interface with eigenvector class
// transform matrix of eigenvectors with gauge array
Eigen::MatrixXcd GaugeField::trafo_ev(const Eigen::MatrixXcd &eig_sys) {
  const ssize_t dim_row = eig_sys.rows();
  const ssize_t dim_col = eig_sys.cols();
  Eigen::MatrixXcd ret(dim_row, dim_col);
  if (omega.shape()[0] == 0)
    build_gauge_array(1);
  // write_gauge_matrices("ev_trafo_log.bin",Omega);

  for (ssize_t nev = 0; nev < dim_col; ++nev) {
    for (ssize_t vol = 0; vol < dim_row; ++vol) {
      int ind_c = vol % 3;
      Eigen::Vector3cd tmp =
          omega[0][ind_c].adjoint() * (eig_sys.col(nev)).segment(ind_c, 3);
      (ret.col(nev)).segment(ind_c, 3) = tmp;
    }  // end loop nev
  }    // end loop vol
  return ret;
}

// Calculate Plaquette for given timeslice

// Plaquette of timeslice
// calculate plaquettes at given point i
double GaugeField::plaque_pnt(const ssize_t mu,
                              const ssize_t nu,
                              const ssize_t vol,
                              const ssize_t t) {
  Eigen::Matrix3cd plaque;
  ssize_t j, k, l;
  double P;
  j = iup[vol][mu];  // convention(up): 0:= n0+1, 1:= n1+1, 2:= n2+1, 3:= n3+1
  k = iup[j][nu];    // convention(down): 0:= n0-1, 1:= n1-1, 2:= n2-1, 3:= n3-1
  l = idown[k][mu];
  Eigen::Matrix3cd a, b, c, d;
  a = tslices.at(t)[vol][mu];  // convention links: 0:= n0+1, 1:= n1+1, 2:= n2+1, 3:= n3+1
  b = tslices.at(t)[j][nu];
  c = tslices.at(t)[l][mu].adjoint();
  d = tslices.at(t)[vol][nu].adjoint();
  plaque = a * (b * (c * d));
  P = 1. / 3. * real(plaque.trace());
  return P;
}

double GaugeField::plaque_ts(const ssize_t t) {
  const ssize_t V3 = Lx * Ly * Lz;

  double plaquette = 0;
  double cnt = 0;
  for (ssize_t vol = 0; vol < V3; ++vol) {
    for (ssize_t dir_1 = 0; dir_1 < 3; ++dir_1) {
      for (ssize_t dir_2 = 0; dir_2 < 3; ++dir_2) {
        // make sure not to run in same direction
        if (dir_1 != dir_2) {
          plaquette += plaque_pnt(dir_1, dir_2, vol, t);
          ++cnt;
        }
      }  // loop dir_2
    }    // loop_dir1
  }      // loop V3
  return (plaquette / cnt);
}

/** Read in gauge field to vector of timeslices */
void GaugeField::read_gauge_field(const ssize_t config_i,
                                  const ssize_t slice_i,
                                  const ssize_t slice_f) {
  const std::string name = config_path + "/conf";
  const auto filename = (boost::format("%s.%04lu") % name % config_i).str();

  std::vector<double> configuration(V_for_lime);

  read_lime_gauge_field_doubleprec_timeslices(
      configuration.data(), filename.c_str(), slice_i, slice_f);
  for (auto t = slice_i; t <= slice_f; ++t) {
    double *const timeslice = &configuration[V_TS * t];
    map_timeslice_to_eigen(t, timeslice);
  }
  std::cout << slice_f + 1 - slice_i << " timeslice(s) read in from config " << config_i
            << std::endl;
}

// Read in the gaugefield, deprecated tmLqcd routine
void GaugeField::read_lime_gauge_field_doubleprec_timeslices(double *gaugefield,
                                                             const char *filename,
                                                             const ssize_t slice_i,
                                                             const ssize_t slice_f) {
  try {
    // const int slice_i = 0;
    // const int slice_f = Lt+1;

    FILE *ifs;
    ssize_t t, x, y, z;
    int status;
    n_uint64_t bytes;
    char *header_type;
    LimeReader *limereader;
    double tmp[72], tmp2[72];
    int words_bigendian;

    printf("reading gauge fields from files:\n");
    //    if(verbose){
    //      printf("reading gauge fields from files:\n");
    //    }
    //    else{
    //      printf("\treading gauge fields:\n");
    //    }

    words_bigendian = big_endian();
    ifs = fopen(filename, "r");
    if (ifs == (FILE *)NULL) {
      fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
      exit(500);
    }
    limereader = limeCreateReader(ifs);
    if (limereader == (LimeReader *)NULL) {
      throw std::runtime_error("Unable to open LimeReader!");
    }
    while ((status = limeReaderNextRecord(limereader)) != LIME_EOF) {
      if (status != LIME_SUCCESS) {
        fprintf(
            stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
        status = LIME_EOF;
        break;
      }
      header_type = limeReaderType(limereader);
      if (strcmp("ildg-binary-data", header_type) == 0)
        break;
    }
    if (status == LIME_EOF) {
      fprintf(stderr, "no ildg-binary-data record found in file %s\n", filename);
      limeDestroyReader(limereader);
      fclose(ifs);
      exit(-2);
    }
    bytes = limeReaderBytes(limereader);
    if (bytes != (n_uint64_t)Lx * Ly * Lz * Lt * 72 * (n_uint64_t)sizeof(double)) {
      if (bytes != (n_uint64_t)Lx * Ly * Lz * Lt * 72 * (n_uint64_t)sizeof(double) / 2) {
        fprintf(stderr,
                "Probably wrong lattice size or precision (bytes=%lu) in file %s "
                "expected %lu\n",
                (n_uint64_t)bytes,
                filename,
                (n_uint64_t)Lx * Ly * Lz * Lt * 72 * (n_uint64_t)sizeof(double));
        fprintf(stderr, "Aborting...!\n");
        fflush(stdout);
        exit(501);
      } else {
        fclose(ifs);
        throw std::runtime_error("Single precision read!");
      }
    }

    bytes = (n_uint64_t)72 * sizeof(double);

    for (t = 0; t < Lt; t++) {
      for (z = 0; z < Lz; z++) {
        for (y = 0; y < Ly; y++) {
          for (x = 0; x < Lx; x++) {
            // check for endianess and reading in data
            // the pointer limereader is internally increased by bytes
            // in the limeReaderReadData function
            if (!words_bigendian) {
              status = limeReaderReadData(tmp, &bytes, limereader);
              byte_swap_assign(tmp2, tmp, 72);
            } else
              status = limeReaderReadData(tmp2, &bytes, limereader);
            // check if reading was successfull
            if (status < 0 && status != LIME_EOR) {
              fprintf(stderr,
                      "LIME read error occured with status = %d while reading file %s!\n "
                      "Aborting...\n",
                      status,
                      filename);
              exit(500);
            }

            // we just want to read in data in the specific range of timeslices
            // must be here because the limereader pointer must be increased
            // correctly
            // could be done with much more performance but might be tricky to
            // do correctly
            if (t < slice_i || t > slice_f)
              continue;

            // copy of link variables from tmp2 into config
            // ILDG has mu-order: x,y,z,t so it is changed here to: t,x,y,z !
            const ssize_t p =
                (ssize_t)(((t - slice_i) * Lx * Lz * Lz + x * Ly * Lz + y * Lz + z) *
                          72);  // position in config
            ssize_t k = 0;
            for (ssize_t mu = 1; mu <= 4; mu++) {  // mu=4 is for the shift of U_t
              ssize_t index;
              if (mu != 4)
                index = p + mu * 18;  // for U_x, U_y and U_z
              else
                index = p;  // U_t is copied into the beginning of
              // the (config+p) array

              for (ssize_t i = 0; i < 3; i++) {
                for (ssize_t j = 0; j < 3; j++) {
                  gaugefield[index + 6 * i + 2 * j] = tmp2[2 * k];
                  gaugefield[index + 6 * i + 2 * j + 1] = tmp2[2 * k + 1];
                  k++;
                }
              }

            }  // loop over mu ends here

          }  // loop over position space ends here
        }
      }
    }
    if (status < 0 && status != LIME_EOR) {
      fprintf(stderr,
              "LIME read error occured with status = %d while reading file %s!\n "
              "Aborting...\n",
              status,
              filename);
      exit(500);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return;

  } catch (std::exception &e) {
    std::cout << e.what()
              << "in: ReadWrite::read_lime_gauge_field_doubleprec_timeslices\n";
    throw e;
  }
}
