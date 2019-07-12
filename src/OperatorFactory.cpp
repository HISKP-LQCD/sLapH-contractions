#include "OperatorsForMesons.hpp"
#include "StopWatch.hpp"
#include "global_data.hpp"

#include <boost/format.hpp>

#include <Eigen/Dense>
#include <iomanip>
#include <cmath>
#include <sstream>

using namespace Eigen;

static ssize_t map_char_to_dir(const char dir) {
  ssize_t integer_dir;
  if (dir == 'x')
    integer_dir = 0;
  if (dir == 'y')
    integer_dir = 1;
  if (dir == 'z')
    integer_dir = 2;
  return integer_dir;
}

/** Displace Matrix of eigenvectors by looping over displacement
 vectors

 v is the Eigensystem at a specific timeslice
 t is the timeslice to take into account
 disp is a vector of displacement pairs
 verbose is an integer controlling the debug level
 The eigenvectors of timeslice t are displaced using the vector of
 displacement pairs. The First entry of each pair specifies the direction
 of the right acting derivative, ">" meaning forward and "<" meaning backward

 @todo Swap order of loops?
 */
Eigen::MatrixXcd displace_eigenvectors(Eigen::MatrixXcd const &v,
                                       GaugeField const &gauge,
                                       ssize_t const t,
                                       DisplacementDirection const disp,
                                       ssize_t const verbose) {
  auto res = v;
  // iterate over displacement vector
  for (const auto &d : disp) {
    // Loop over all eigenvectors in v
    for (int ev = 0; ev < v.cols(); ++ev) {
      // Multiply eigenvector with according gauge matrix
      for (int spatial_ind = 0; spatial_ind < v.rows() / 3; ++spatial_ind) {
        res.col(ev).segment(3 * spatial_ind, 3) =
            (d.first == '>')
                ? gauge.forward_uv(res.col(ev), t, spatial_ind, map_char_to_dir(d.second))
                : gauge.backward_uv(
                      res.col(ev), t, spatial_ind, map_char_to_dir(d.second));
      }  // End spatial loop

    }  // end eigenvector loop
  }
  return res;
}

/** Creates a two-dimensional vector containing the momenta for the operators
 *
 *  @param[in] Lx, Ly, Lz      Lattice extent in spatial directions
 *  @param[in] vdaggerv_lookup Contains the momenta
 *  @param[in,out] momentum    Two dimensional array where the momenta are
 *                             stored
 */
void create_momenta(ssize_t const Lx,
                    ssize_t const Ly,
                    ssize_t const Lz,
                    std::vector<VdaggerVQuantumNumbers> const &vdaggerv_lookup,
                    array_cd_d2 &momentum) {
  static const std::complex<double> I(0.0, 1.0);

  /** To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
   *  quantum number id in op_VdaggerV will be used. The rest can be obtained
   *  by adjoining
   */
  for (const auto &op : vdaggerv_lookup) {
    // op_VdaggerV contains the index of one (redundancy) op_Corr which
    // allows to deduce the quantum numbers (momentum)
    const double ipx = op.momentum[0] * 2. * M_PI / static_cast<double>(Lx);
    const double ipy = op.momentum[1] * 2. * M_PI / static_cast<double>(Ly);
    const double ipz = op.momentum[2] * 2. * M_PI / static_cast<double>(Lz);
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

std::vector<Complex> roots_of_unity(int momentum, ssize_t size) {
  Complex const I(0.0, 1.0);

  double const two_pi_p_L = momentum * 2.0 * M_PI / static_cast<double>(size);

  std::vector<Complex> roots(size);
  for (ssize_t i = 0; i < size; ++i) {
    roots[i] = exp(-I * two_pi_p_L * static_cast<double>(i));
  }

  return roots;
}

void create_momenta_new(ssize_t const Lx,
                        ssize_t const Ly,
                        ssize_t const Lz,
                        std::vector<VdaggerVQuantumNumbers> const &vdaggerv_lookup,
                        array_cd_d2 &momentum) {
  /** To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
   *  quantum number id in op_VdaggerV will be used. The rest can be obtained
   *  by adjoining
   */
  for (auto const &op : vdaggerv_lookup) {
    auto const phases_x = roots_of_unity(op.momentum[0], Lx);
    auto const phases_y = roots_of_unity(op.momentum[1], Ly);
    auto const phases_z = roots_of_unity(op.momentum[2], Lz);

    // Calculate $\vec{p} \cdot \vec{x}$ for all $\vec{x}$ on the lattice.
    for (ssize_t x = 0; x < Lx; ++x) {
      auto const xH = x * Ly * Lz;
      auto const phase_x = phases_x[x];

      for (ssize_t y = 0; y < Ly; ++y) {
        auto const xHyH = xH + y * Lz;
        auto const phase_xy = phase_x * phases_y[y];

        for (ssize_t z = 0; z < Lz; ++z) {
          auto const phase_xyz = phase_xy * phases_z[z];
          momentum[op.id][xHyH + z] = phase_xyz;
        }
      }
    }
  }
}

namespace {

void write_vdaggerv(const std::string &pathname,
                    const std::string &filename,
                    const Eigen::MatrixXcd &Vt) {
  // writing the data
  std::ofstream file((pathname + filename).c_str(),
                     std::ofstream::binary | std::ofstream::trunc);

  if (file.is_open()) {
    std::cout << "\twriting VdaggerV to file:" << pathname + filename << std::endl;
    // buffer for writing
    std::vector<Complex> eigen_vec(Vt.size());
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

}  // namespace

/**
 * @param Lt, Lx, Ly, Lz  Temporal and spatial lattice extent
 * @param nb_ev           Number of eigenvectors
 * @param dilE            Number of diluted blocks in eigenvector space
 * @param operator_lookuptable ?
 * @param handling_vdaggerv
 * @param path_vdaggerv
 *
 * The initialization of the container attributes of OperatorsForMesons
 * is done in the member initializer list of the constructor. The allocation
 * of heap memory is delegated to boost::multi_array::resize
 */
OperatorFactory::OperatorFactory(const ssize_t Lt,
                                 const ssize_t Lx,
                                 const ssize_t Ly,
                                 const ssize_t Lz,
                                 const ssize_t nb_ev,
                                 const ssize_t dilE,
                                 const OperatorLookup &operator_lookuptable,
                                 const std::string &handling_vdaggerv,
                                 const std::string &path_vdaggerv,
                                 const std::string &path_config,
                                 const HypPars &hyp_parameters)
    : vdaggerv(),
      momentum(),
      operator_lookuptable(operator_lookuptable),
      Lt(Lt),
      Lx(Lx),
      Ly(Ly),
      Lz(Lz),
      nb_ev(nb_ev),
      dilE(dilE),
      handling_vdaggerv(handling_vdaggerv),
      path_vdaggerv(path_vdaggerv),
      path_config(path_config),
      hyp_parameters(hyp_parameters) {
  // resizing containers to their correct size
  vdaggerv.resize(boost::extents[operator_lookuptable.vdaggerv_lookup.size()][Lt]);

  // the momenta only need to be calculated for a subset of quantum numbers
  // (see VdaggerV::build_vdaggerv)
  momentum.resize(
      boost::extents[operator_lookuptable.vdaggerv_lookup.size()][Lx * Ly * Lz]);
  create_momenta(Lx, Ly, Lz, operator_lookuptable.vdaggerv_lookup, momentum);

  std::cout << "\tMeson operators initialised" << std::endl;
}

static inline void kernel_compute_vdaggerv(const ssize_t dim_row,
                                           const ssize_t nb_ev,
                                           const int config, 
                                           const int phase_t_start_idx,
                                           const int phase_nb_t_slices,
                                           const EigenVector &V_t, 
                                           const array_cd_d2 &momentum, 
                                           const GlobalData & gd, 
                                           const OperatorLookup &operator_lookuptable, 
                                           array_Xcd_d2_eigen & vdaggerv, 
                                           const GaugeField &gauge ){

  Eigen::MatrixXcd W_t;
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row/3);

  const int id_unity = operator_lookuptable.index_of_unity;
  
  const int old_eigen_threads = Eigen::nbThreads();
  Eigen::setNbThreads(gd.nb_vdaggerv_eigen_threads);

  StopWatch vdaggerv_check_watch("VdaggerV check (note: nb_vdagger_eigen_threads used)", 1);
  vdaggerv_check_watch.start();
  for(ssize_t i = 0; i < phase_nb_t_slices; ++i){
    bool test_fail = V_t.test_trace_sum(i, false);
    if(test_fail){
      std::stringstream message;
      message << "Eigenvector verification failed at config: " << config << " time slice: " <<
         phase_t_start_idx + i << std::endl;
      throw std::runtime_error( message.str() );
    }
  }
  vdaggerv_check_watch.stop();
  vdaggerv_check_watch.print();

  //Look up table for the off-diagonal part
  //in order to efficiently parallelize te
  //loop over the off diagonal entries of vdaggerv
  std::vector<int> looplookuptable;
  for (ssize_t ncol1=0; ncol1<(V_t[0]).cols(); ++ncol1) {
    for (ssize_t ncol2=ncol1+1; ncol2<(V_t[0]).cols(); ++ncol2) {
      looplookuptable.push_back(ncol2+(V_t[0]).cols()*ncol1);
    }
  }


  for(ssize_t i = 0; i < phase_nb_t_slices; ++i){
    const ssize_t t = phase_t_start_idx + i;
    for(const auto &op : operator_lookuptable.vdaggerv_lookup){
      if( op.id != id_unity ){
       if(!op.displacement.empty()){
         W_t.noalias() = displace_eigenvectors(V_t[i], gauge, t, op.displacement, 1);
         vdaggerv[op.id][t] = V_t[i].adjoint() * W_t;
       } 
      }
      if ((op.id == id_unity) && (op.displacement.empty())){
       vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
      }
    }
    //Computation of Vdaggerv: the main idea is 
    //     (a) reduce the size of dense matrix multiplication by a factor of three
    //     (b) parallelize not in the x spatial coordinate, but in the index
    //         of vdaggerv itself
    //     (c) taking advantage of the fact that vdaggerv is a symmetric matrix
    //Computing the diagonal part of VdaggerV
    # pragma omp parallel for num_threads(gd.nb_vdaggerv_eigen_threads) schedule(static)
    for(ssize_t ncol = 0; ncol < V_t[i].cols() ; ++ncol){
      //Trick: reshape the eigenvectors in order to take the scalar products in color space
      MatrixXcd vnew=V_t[i];//I would like to avoid this copy, however
      //Eigen::Map<Eigen::MatrixXcd> Reshaped((V_t)[i].col(ncol).data(), 3, (V_t[i].col(ncol)).size()/3);
      //does not work
      Eigen::Map<Eigen::MatrixXcd> Reshaped(vnew.col(ncol).data(), 3, (V_t[i].col(ncol)).size()/3);
      //Taking the scalar products
      MatrixXcd xspaceresults=Reshaped.colwise().squaredNorm();
      //Now we do the multiplication for each possible momentum
      for(const auto &op : operator_lookuptable.vdaggerv_lookup){
        if( op.id != id_unity ){
         if(op.displacement.empty()){
          for (ssize_t x = 0; x < dim_row/3; ++x) {
            mom(x) = momentum[op.id][x];
          }
          vdaggerv[op.id][t](ncol,ncol)=xspaceresults.row(0).dot(mom);
         }
        }
      }
    }
    //Computing the off-diagonal part of VdaggerV
    # pragma omp parallel for num_threads(gd.nb_vdaggerv_eigen_threads) schedule(static)
    for (ssize_t nindex=0; nindex<(V_t[i].cols()*V_t[i].cols()-V_t[i].cols())/2;++nindex){
      ssize_t ncol=looplookuptable[nindex]/V_t[i].cols();
      ssize_t nrow=looplookuptable[nindex]%V_t[i].cols();
      //We apply the same trick as for the diagonal elements: reshape
      //the eigenvector and then compute elementwise product of the two
      //matrix
      MatrixXcd vnew=V_t[i];
      Eigen::Map<Eigen::MatrixXcd> Reshaped1(vnew.col(ncol).data(), 3, (V_t[i].col(ncol)).size()/3);
      Eigen::Map<Eigen::MatrixXcd> Reshaped2(vnew.col(nrow).data(), 3, (V_t[i].col(nrow)).size()/3);
      //Taking the elementwise product
      MatrixXcd xspaceresults = (Reshaped2.array().conjugate()*Reshaped1.array());
      //MatrixXcd xspaceresults = Reshaped2.cwiseProduct(Reshaped1);
      //Here we do not take the SquaredNorm, because the multiplication was
      //already done 
      for(const auto &op : operator_lookuptable.vdaggerv_lookup){
        if( op.id != id_unity ){
         if(op.displacement.empty()){
          for (ssize_t x = 0; x < dim_row/3; ++x) {
            mom(x) = momentum[op.id][x];
          }
          vdaggerv[op.id][t](ncol,nrow)= xspaceresults.colwise().sum().dot(mom);
          vdaggerv[op.id][t](nrow,ncol)= vdaggerv[op.id][t](ncol, nrow);
         }
        }
      }
    }
  }
  Eigen::setNbThreads(old_eigen_threads);
}

void OperatorFactory::build_vdaggerv(const std::string &filename, const int config, const GlobalData & gd) {
  const ssize_t dim_row = 3 * Lx * Ly * Lz;
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for writing
  std::string const full_path =
      (boost::format("/%s/cnfg%04d/") % path_vdaggerv % config).str();

  // check if directory exists
  if (( (handling_vdaggerv == "write") || ( handling_vdaggerv == "only_vdaggerv_compute_save" ) )&& access(full_path.c_str(), 0) != 0) {
    std::cout << "\tdirectory " << full_path.c_str()
              << " does not exist and will be created";
    boost::filesystem::path dir(full_path.c_str());
    if (!boost::filesystem::create_directories(dir))
      std::cout << "\tSuccess" << std::endl;
    else
      std::cout << "\tFailure" << std::endl;
  }

  StopWatch swatch("Eigenvector and Gauge I/O and VdaggerV construction (note: various thread numbers used)",
                   1);
  swatch.start();
  // resizing each matrix in vdaggerv
  // TODO: check if it is better to use for_each and resize instead of std::fill
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

  // Read gauge field only if it is needed.
  /** @todo might be useful to parallelize */
  std::unique_ptr<GaugeField> gauge(nullptr);
  if (operator_lookuptable.need_gaugefield) {
    // If parameters for smearing are set, smear operator
    gauge.reset(new GaugeField(Lt, Lx, Ly, Lz, path_config, 0, Lt - 1, 4));
    gauge->read_gauge_field(config, 0, Lt - 1);
    if (hyp_parameters.iterations > 0) {
      const double alpha1 = hyp_parameters.alpha1;
      const double alpha2 = hyp_parameters.alpha2;
      const size_t iter = hyp_parameters.iterations;
      for (ssize_t t = 0; t < Lt; ++t) {
        gauge->smearing_hyp(t, alpha1, alpha2, iter);
      }
    }
  }

  // memory for eigensystems on 'nb_evec_read_threads' time slices
  StopWatch malloc_watch("V_t memory allocation", 1);
  malloc_watch.start();
  EigenVector V_t(gd.nb_evec_read_threads, dim_row, nb_ev);
  malloc_watch.stop();
  malloc_watch.print();

  // how many outer iterations do we need to work off all time slices in a stepping 
  // of 'nb_evec_read_threads'
  const ssize_t nphases = (ssize_t)ceil((double)Lt / gd.nb_evec_read_threads);
 
  // starting time slice 
  ssize_t phase_t_start_idx = 0;
 
  // loop over the phases
  StopWatch phase_watch("Outer loop phase", 1);
  for(ssize_t iphase = 0; iphase < nphases; iphase++){
    phase_watch.start();
    std::cout << "Phase " << iphase+1 << " of " << nphases << std::endl;

    ssize_t t_end = phase_t_start_idx + gd.nb_evec_read_threads;
    ssize_t phase_nb_t_slices = gd.nb_evec_read_threads;
    if( t_end >= Lt ){
      phase_nb_t_slices = Lt - phase_t_start_idx;
      t_end = Lt;
    }

    StopWatch ev_io_watch("Threaded Eigenvector I/O", phase_nb_t_slices);
    // perform thread-parallel I/O with the thread number chosen so as to not exceed node
    // memory limits
    #pragma omp parallel num_threads(phase_nb_t_slices)
    {
      ev_io_watch.start();
      #pragma omp for schedule(static)
      for(ssize_t i = 0; i < phase_nb_t_slices; ++i){
        ssize_t t = phase_t_start_idx + i;
        auto const inter_name = (boost::format("%s%03d") % filename % t).str();
        // we call with verbose=0 because we want to run verification in the more efficient mode
        // down below, when 'nb_vdaggerv_eigen_threads' are in use to perform VdV.trace and VdV.sum
        V_t.read_eigen_vector(inter_name.c_str(), i, 0, false);

        // for small lattices, it is most efficient to trivially parallelize the
        // computation of vdaggerv, running with 'nb_vdaggerv_eigen_threads == 1'
        // and as many reading threads as is efficient
        if (gd.nb_vdaggerv_eigen_threads == 1) {
          kernel_compute_vdaggerv(dim_row,
                                  nb_ev,
                                  config,
                                  phase_t_start_idx,
                                  phase_nb_t_slices,
                                  V_t,
                                  momentum,
                                  gd,
                                  operator_lookuptable,
                                  vdaggerv,
                                  *gauge);
        }

      }
      ev_io_watch.stop();
    }
    ev_io_watch.print();

    if (gd.nb_vdaggerv_eigen_threads != 1) {
      // for larget lattices, one can avoid exceeding memory limits by doing threaded I/O above
      // with 'nb_evec_read_threads' while doing the dense linear algebra with 'nb_vdaggerv_eigen_threads'
      // down here, probably reasonably efficiently.
      // Unfortunately, Eigen self-limits the level of parallelism in the VdaggerV computation
      // For instance, on a 32c64 lattice with 220 eigenvectors, it will limit itself to six threads
      // no matter what has been set of nb_vdaggerv_eigen_threads
      // On large lattices, however, allowing all cores to be used leads to lower total time
      // for build_vdaggerv
      kernel_compute_vdaggerv(dim_row,
                              nb_ev,
                              config,
                              phase_t_start_idx,
                              phase_nb_t_slices,
                              V_t,
                              momentum,
                              gd,
                              operator_lookuptable,
                              vdaggerv,
                              *gauge );
    }
    
    phase_t_start_idx += phase_nb_t_slices;

    phase_watch.stop();
    phase_watch.print();
  } // for(iphase)
  

  // perform thread-parallel I/O to write out VdaggerV
  if( ( handling_vdaggerv == "write") || ( handling_vdaggerv == "only_vdaggerv_compute_save" )){
    StopWatch vdaggerv_io_watch("VdaggerV thread-parallel writing",
                                gd.nb_evec_read_threads);
    #pragma omp parallel num_threads(gd.nb_evec_read_threads)
    {
      vdaggerv_io_watch.start();
      #pragma omp for schedule(dynamic)
      for(ssize_t t = 0; t < Lt; ++t){
        for(const auto &op : operator_lookuptable.vdaggerv_lookup){
          if( op.id != id_unity ){
            std::string momentum_string = std::to_string(op.momentum[0]) +
                                          std::to_string(op.momentum[1]) +
                                          std::to_string(op.momentum[2]);
            std::string displacement_string = to_string(op.displacement);
            std::string outfile =
                (boost::format("operators.%04d.p_%s.d_%s.t_%03d") % config %
                 momentum_string % displacement_string % (int)t)
                    .str();
            write_vdaggerv(full_path, std::string(outfile), vdaggerv[op.id][t]);
          }
        }
      }
      vdaggerv_io_watch.stop();
    } // pragma omp parallel
    vdaggerv_io_watch.print();
  }

  swatch.stop();
  swatch.print();
  is_vdaggerv_set = true;
}

void OperatorFactory::read_vdaggerv(const int config) {
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for reading
  auto const full_path =
      (boost::format("/%s/cnfg%04d/operators.%04d") % path_vdaggerv % config % config)
          .str();

  // resizing each matrix in vdaggerv
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));
  StopWatch swatch("VdaggerV I/O");

#pragma omp parallel
  {
    swatch.start();
#pragma omp for schedule(dynamic)
    for (ssize_t t = 0; t < Lt; ++t) {
      for (const auto &op : operator_lookuptable.vdaggerv_lookup) {
        // For zero momentum and displacement VdaggerV is the unit matrix, thus
        // the calculation is not performed
        if (op.id != id_unity) {
          // creating full filename for vdaggerv and reading them in
          std::string dummy = full_path + ".p_" + std::to_string(op.momentum[0]) +
                              std::to_string(op.momentum[1]) +
                              std::to_string(op.momentum[2]);

          auto const infile = (boost::format("%s_.t_%03d") % dummy % t).str();

          // writing the data
          std::ifstream file(infile, std::ifstream::binary);

          if (file.is_open()) {
            std::cout << "\treading VdaggerV from file:" << infile << std::endl;

            // buffer for reading
            std::vector<Complex> eigen_vec(vdaggerv[op.id][t].size());
            file.read(reinterpret_cast<char *>(&eigen_vec[0]),
                      vdaggerv[op.id][t].size() * sizeof(Complex));
            for (ssize_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (ssize_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(ncol * vdaggerv[op.id][t].rows() + nrow);
              }
            }
            if (!file.good()) {
              std::ostringstream oss;
              oss << "Problems while reading from " << infile;
              std::runtime_error(oss.str());
            }
            file.close();
          } else {
            std::ostringstream oss;
            oss << "Can't open " << infile;
            std::runtime_error(oss.str());
          }
        } else  // zero momentum
          vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
      }
    }  // loop over time
    swatch.stop();
  }  // pragma omp parallel ends here

  swatch.print();
  is_vdaggerv_set = true;
}

void OperatorFactory::read_vdaggerv_liuming(const int config) {
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for reading
  auto const full_path = (boost::format("/%s/VdaggerV.") % path_vdaggerv).str();

  // resizing each matrix in vdaggerv
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

  StopWatch swatch("Liuming VdaggerV I/O");
#pragma omp parallel
  {
    swatch.start();
    //  #pragma omp for schedule(dynamic)
    //    for(const auto& op : operator_lookuptable.vdaggerv_lookup){
#pragma omp for schedule(dynamic)
    for (ssize_t i = 0; i < ssize(operator_lookuptable.vdaggerv_lookup); ++i) {
      const auto op = (operator_lookuptable.vdaggerv_lookup[i]);
      // For zero momentum and displacement VdaggerV is the unit matrix, thus
      // the calculation is not performed
      if (op.id != id_unity) {
        // creating full filename for vdaggerv and reading them in
        // both possibilities must be checked
        std::string dummy1 = full_path + "p" + std::to_string(-op.momentum[0]) + "p" +
                             std::to_string(-op.momentum[1]) + "p" +
                             std::to_string(-op.momentum[2]) + ".conf";
        auto const infile1 = (boost::format("%s%04d") % dummy1 % config).str();
        std::ifstream file1(infile1, std::ifstream::binary);

        // second possibility for a name
        std::string dummy2 = full_path + "p" + std::to_string(op.momentum[0]) + "p" +
                             std::to_string(op.momentum[1]) + "p" +
                             std::to_string(op.momentum[2]) + ".conf";
        auto const infile2 = (boost::format("%s%04d") % dummy2 % config).str();
        std::ifstream file2(infile2, std::ifstream::binary);

        if (file1.is_open()) {
          std::cout << "\treading VdaggerV from file:" << infile1 << std::endl;
          for (ssize_t t = 0; t < Lt; ++t) {
            // buffer for reading
            std::vector<Complex> eigen_vec(vdaggerv[op.id][t].size());
            file1.read(reinterpret_cast<char *>(&eigen_vec[0]),
                       vdaggerv[op.id][t].size() * sizeof(Complex));
            for (ssize_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (ssize_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(nrow * vdaggerv[op.id][t].cols() + ncol);
              }
            }
            vdaggerv[op.id][t].adjointInPlace();
            if (!file1.good()) {
              std::ostringstream oss;
              oss << "Problems while reading from " << infile1;
              std::runtime_error(oss.str());
            }
          }  // loop over time
          file1.close();
        } else if (file2.is_open()) {
          std::cout << "\treading VdaggerV from file:" << infile2 << std::endl;
          for (ssize_t t = 0; t < Lt; ++t) {
            // buffer for reading
            std::vector<Complex> eigen_vec(vdaggerv[op.id][t].size());
            file2.read(reinterpret_cast<char *>(&eigen_vec[0]),
                       vdaggerv[op.id][t].size() * sizeof(Complex));
            for (ssize_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (ssize_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(nrow * vdaggerv[op.id][t].cols() + ncol);
              }
            }
            //            // @todo Check whether that must be adjoint before file 1 is
            //            closed
            //            // (master branch)
            //            vdaggerv[op.id][t].adjointInPlace();
            if (!file2.good()) {
              std::ostringstream oss;
              oss << "Problems while reading from " << infile2;
              std::runtime_error(oss.str());
            }
          }  // loop over time
          file2.close();
        } else {
          std::ostringstream oss;
          oss << "can't open " << infile1 << " NOR " << infile2;
          std::runtime_error(oss.str());
        }
      } else  // zero momentum
        for (ssize_t t = 0; t < Lt; ++t)
          vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
    }
    swatch.stop();
  }  // pragma omp parallel ends here

  swatch.print();
  is_vdaggerv_set = true;
}

/**
 *  @param filename The name to write to / read from the V^\dagger V operators
 *  @param rnd_vec  The random vector
 *  @param config   The configuration number to be read. Unused if
 *                  handling_vdaggerv is "build" or "write"
 *
 *  Behavior of this function depends on handling_vdaggerv flag.
 *  - "only_vdaggerv_compute_save" Only the vdaggerv objects are construncted and written out
 *  - "read" | "liuming" The operators are read in the corresponding format.
 *  - "build"            The operators are constructed from the eigenvectors
 *  - "write"            The operators are constructed and additionaly written
 *                       out.
 */
void OperatorFactory::create_operators(const std::string &filename,
                                       const RandomVector &rnd_vec,
                                       const int config,
                                       const GlobalData & gd) {
  is_vdaggerv_set = false;
  if (handling_vdaggerv == "write" || handling_vdaggerv == "build" || handling_vdaggerv == "only_vdaggerv_compute_save")
    build_vdaggerv(filename, config, gd);
  else if (handling_vdaggerv == "read")
    read_vdaggerv(config);
  else if (handling_vdaggerv == "liuming")
    read_vdaggerv_liuming(config);
  else {
    throw std::runtime_error("The flag handling_vdaggerv in input file is wrong!");
  }
}

/**
 *  E.g. after building Quarkline Q2, vdaggerv is no longer needed and can be
 *  deleted to free up space
 *
 *  Resizes vdaggerv to 0
 */
void OperatorFactory::free_memory_vdaggerv() {
  std::for_each(vdaggerv.origin(),
                vdaggerv.origin() + vdaggerv.num_elements(),
                [](Eigen::MatrixXcd m) { m.resize(0, 0); });
}
