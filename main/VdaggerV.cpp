#include "Correlators.hpp"
#include "git.hpp"
#include "global_data.hpp"
#include "global_data_build_IO_names.hpp"
#include "timings.hpp"

#include <omp.h>

#include <iostream>
#include "OperatorsForMesons.hpp"

/*
#include "EigenVector.hpp"
#include "typedefs.hpp"

#include <omp.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <iomanip>
#include <iostream>

void read_eigen_vector(const std::string &filename,
                       const ssize_t verbose,
                       Eigen::MatrixXcd &Vt) {
  // buffer for read in
  vec eigen_vec(Vt.rows());
  //  std::cout << "\tReading eigenvectors from files:" << filename << std::endl;

  // setting up file
  std::ifstream infile(filename, std::ifstream::binary);
  if (infile) {
    for (ssize_t ncol = 0; ncol < Vt.cols(); ++ncol) {
      infile.read((char *)&(eigen_vec[0]), 2 * Vt.rows() * sizeof(double));
      for (ssize_t nrow = 0; nrow < Vt.rows(); ++nrow) {
        (Vt)(nrow, ncol) = eigen_vec[nrow];
      }
    }
  } else {
    throw std::runtime_error("Eigenvector file does not exist!");
  }
  infile.close();

  // small test of trace and sum over the eigen vector matrix!
  if (verbose) {
    std::cout << "trace of V^d*V"
              << ":\t" << (Vt.adjoint() * Vt).trace() << std::endl;
    std::cout << "sum over all entries of V^d*V"
              << ":\t" << (Vt.adjoint() * Vt).sum() << std::endl;
  }
}
void write_vdaggerv(const std::string &pathname,
                    const std::string &filename,
                    const Eigen::MatrixXcd &Vt) {
  // check if directory exists
  if (access(pathname.c_str(), 0) != 0) {
    std::cout << "\tdirectory " << pathname.c_str()
              << " does not exist and will be created";
    boost::filesystem::path dir(pathname.c_str());
    if (!boost::filesystem::create_directories(dir))
      std::cout << "\tSuccess" << std::endl;
    else
      std::cout << "\tFailure" << std::endl;
  }
  // writing the data
  std::ofstream file((pathname + filename).c_str(),
                     std::ofstream::binary | std::ofstream::trunc);

  if (file.is_open()) {
    std::cout << "\twriting VdaggerV to file:" << pathname + filename << std::endl;
    // buffer for writing
    vec eigen_vec(Vt.size());
    for (ssize_t ncol = 0; ncol < Vt.cols(); ncol++) {
      for (ssize_t nrow = 0; nrow < Vt.rows(); nrow++) {
        eigen_vec.at(ncol * Vt.rows() + nrow) = (Vt)(nrow, ncol);
      }
    }
    file.write(reinterpret_cast<const char *>(&eigen_vec[0]),
               Vt.size() * sizeof(Complex));
    if (!file.good())
      std::cout << "Problems while write to " << (pathname + filename).c_str()
                << std::endl;
    file.close();
  } else
    std::cout << "can't open " << (pathname + filename).c_str() << std::endl;
}

// creates a two-dimensional vector containing the momenta for the operators
// input: Lx, Ly, Lz      -> lattice extend in x, y, and z direction
//        vdaggerv_lookup -> contains the momenta
// output: momentum -> two dimensional array where the momenta are stored
void create_momenta(const ssize_t Lx,
                    const ssize_t Ly,
                    const ssize_t Lz,
                    const std::vector<VdaggerVQuantumNumbers> &vdaggerv_lookup,
                    array_cd_d2 &momentum) {
  static const std::complex<double> I(0.0, 1.0);

  // To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
  // quantum number id in op_VdaggerV will be used. The rest can be obtained
  // by adjoining
  for (const auto &op : vdaggerv_lookup) {
    // op_VdaggerV contains the index of one (redundancy) op_Corr which
    // allows to deduce the quantum numbers (momentum)
    const double ipx = op.momentum[0] * 2. * M_PI / (double)Lx;
    const double ipy = op.momentum[1] * 2. * M_PI / (double)Ly;
    const double ipz = op.momentum[2] * 2. * M_PI / (double)Lz;
    // calculate \vec{p} \cdot \vec{x} for all \vec{x} on the lattice
    for (int x = 0; x < Lx; ++x) {
      const int xH = x * Ly * Lz;   // helper variable
      const double ipxH = ipx * x;  // helper variable
      for (int y = 0; y < Ly; ++y) {
        const int xHyH = xH + y * Lz;            // helper variable
        const double ipxHipyH = ipxH + ipy * y;  // helper variable
        for (int z = 0; z < Lz; ++z) {
          // multiply \vec{p} \cdot \vec{x} with complex unit and exponentiate
          momentum[op.id][xHyH + z] = exp(-I * (ipxHipyH + ipz * z));
        }
      }
    }  // loops over spatial vectors end here
  }    // loop over redundant quantum numbers ends here
}

void build_and_write_vdaggerv(const ssize_t Lt,
                              const ssize_t Lx,
                              const ssize_t Ly,
                              const ssize_t Lz,
                              const ssize_t nb_ev,
                              const OperatorLookup &operator_lookuptable,
                              const std::string &filename_ev,
                              const std::string &pathname_vdaggerv,
                              const std::string &filename_vdaggerv) {
  clock_t t2 = clock();
  const ssize_t dim_row = 3 * Lx * Ly * Lz;
  const ssize_t id_unity = operator_lookuptable.index_of_unity;

  array_cd_d2 momentum;
  momentum.resize(
      boost::extents[operator_lookuptable.vdaggerv_lookup.size()][Lx * Ly * Lz]);
  create_momenta(Lx, Ly, Lz, operator_lookuptable.vdaggerv_lookup, momentum);

#pragma omp parallel
  {
    Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);
    Eigen::MatrixXcd vdaggerv = Eigen::MatrixXcd::Zero(nb_ev, nb_ev);
    Eigen::MatrixXcd V_t = Eigen::MatrixXcd::Zero(dim_row, nb_ev);

#pragma omp for schedule(dynamic)
    for (ssize_t t = 0; t < Lt; ++t) {
      // creating full filename for eigenvectors and reading them in
      char infile[200];
      sprintf(infile, "%s.%03d", filename_ev.c_str(), (int)t);
      read_eigen_vector(infile, 0, V_t);  // reading eigenvectors

      // VdaggerV is independent of the gamma structure and momenta connected by
      // sign flip are related by adjoining VdaggerV. Thus the expensive
      // calculation must only be performed for a subset of quantum numbers given
      // in op_VdaggerV.
      for (const auto &op : operator_lookuptable.vdaggerv_lookup) {
        // For zero momentum and displacement VdaggerV is the unit matrix, thus
        // the calculation is not performed
        if (op.id != id_unity) {
          // momentum vector contains exp(-i p x). Divisor 3 for colour index.
          // All three colours on same lattice site get the same momentum.
          for (ssize_t x = 0; x < dim_row; ++x) {
            mom(x) = momentum[op.id][x / 3];
          }
          vdaggerv = V_t.adjoint() * mom.asDiagonal() * V_t;

          // write operators to file
          std::string dummy = filename_vdaggerv + ".p_" + std::to_string(op.momentum[0]) +
                              std::to_string(op.momentum[1]) +
                              std::to_string(op.momentum[2]);
          char outfile[200];
          sprintf(outfile, "%s_.t_%03d", dummy.c_str(), (int)t);
          write_vdaggerv(pathname_vdaggerv, std::string(outfile), vdaggerv);
        }
      }
    }  // loop over time
  }    // pragma omp parallel ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed
            << ((float)t2) / CLOCKS_PER_SEC << " seconds" << std::endl;
}
*/
int main(int ac, char *av[]) {
  // Some variables definitions which should be read from infile!
  //
   std::cout << "This is sLapH-contractions:\n"
            << "  git branch " << git_refspec << "\n"
            << "  git revision " << git_sha1 << "\n"
            << "  git state " << git_changes << "\n"
            << "  compiled by " << git_user << " on " << git_host << "\n"
            << "  running with up to " << omp_get_max_threads() << " OpenMP threads\n";

  // reading global parameters from input file
  GlobalData gd;
  read_parameters(gd, ac, av);
  std::cout<<"Parameters are read"<<std::endl;
    

  // initialization of OMP paralization
  Eigen::initParallel();
  Eigen::setNbThreads(gd.nb_eigen_threads);
  
  // Creating lookuptable for Operator construction by hand. TODO: Include
  // this in input file handling!
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup;
  // emplacing back all the momenta
  //std::array<int, 3> displacement = {0, 0, 0};
  std::vector<std::pair<char, char>> displacement;
  ssize_t id = 0;
  for (int p1 = gd.max_momentum; p1 >= -gd.max_momentum; p1--)
    for (int p2 = gd.max_momentum; p2 >= -gd.max_momentum; p2--)
      for (int p3 = gd.max_momentum; p3 >= -gd.max_momentum; p3--){
        int maxQsq=gd.max_momentum*gd.max_momentum;
        if (p1 * p1 + p2 * p2 + p3 * p3 <= maxQsq && !(p1 == 0 && p2 == 0 && p3 == 0)) {
          std::array<int, 3> momentum = {p1, p2, p3};
          std::cout << p1 << "\t" << p2 << "\t" << p3 << std::endl;
          auto it = std::find_if(
              vdaggerv_lookup.begin(),
              vdaggerv_lookup.end(),
              [&momentum](VdaggerVQuantumNumbers vdv_qn) {
                const std::array<int, 3> pm = {-momentum[0], -momentum[1], -momentum[2]};
                if (vdv_qn.momentum == pm)
                  return true;
                else
                  return false;
              });
          if (!(it != vdaggerv_lookup.end())) {
            vdaggerv_lookup.emplace_back(
                id, momentum, displacement);
            id++;
          }
        }
  }
  OperatorLookup operator_lookuptable;
  operator_lookuptable.vdaggerv_lookup = vdaggerv_lookup;
  operator_lookuptable.index_of_unity = 0;


  // Loop over all configurations stated in the infile -------------------------
  for (ssize_t config_i = gd.start_config; config_i <= gd.end_config;
       config_i += gd.delta_config) {
    std::cout << "\nprocessing configuration: " << config_i << "\n\n";

    build_IO_names(gd, config_i);
    OperatorFactory meson_operators(gd.Lt,
                                  gd.Lx,
                                  gd.Ly,
                                  gd.Lz,
                                  gd.number_of_eigen_vec,
                                  (gd.quarks)[0].number_of_dilution_E,
                                  gd.operator_lookuptable,
                                  gd.handling_vdaggerv,
                                  gd.path_vdaggerv,
                                  gd.path_config,
                                  gd.hyp_parameters);

    meson_operators.build_vdaggerv(gd.filename_eigenvectors, config_i, gd);


  }
}
