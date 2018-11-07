#include "EigenVector.hpp"

#include <boost/format.hpp>

void EigenVector::write_eigen_vector(const std::string &filename,
                                     const ssize_t t,
                                     const ssize_t verbose){
  // setting up file
  std::ofstream outfile(filename, std::ofstream::binary);
  if (outfile) {
      std::streamsize  begin = outfile.tellp();
      std::streamsize eigsys_bytes =2*V[t].rows()*V[t].cols()*sizeof(double); 
      outfile.write(reinterpret_cast<char*> (V[t].data()), eigsys_bytes);
      std::streamsize end = outfile.tellp();
      if ( (end - begin)/eigsys_bytes != 1 ){
        std::ostringstream oss;
        oss << "Timeslice:  " << t << ". Error: write incomplete, exiting. "
            << (end - begin) << " bytes instead of expected " << eigsys_bytes << " bytes";
        throw std::runtime_error(oss.str());
      } 
  } else {
    throw std::runtime_error("Eigenvector file does not exist!");
  }
  outfile.close();

  // small test of trace and sum over the eigen vector matrix!
  if (verbose) {
    std::cout << "trace of V^d*V"
              << ":\t" << (V[t].adjoint() * V[t]).trace() << std::endl;
    std::cout << "sum over all entries of V^d*V"
              << ":\t" << (V[t].adjoint() * V[t]).sum() << std::endl;
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void EigenVector::read_eigen_vector(const std::string &filename,
                                    const ssize_t t,
                                    const ssize_t verbose,
                                    const bool mock) {
  if(!mock){
    // buffer for read in
    std::vector<Complex> eigen_vec(V[t].rows());
    std::cout << "\tReading eigenvectors from files:" << filename << std::endl;

    // setting V[t] to zero
    V[t].setZero();

    //// setting up file
    std::ifstream infile(filename, std::ifstream::binary);
    if (infile) {
      for (ssize_t ncol = 0; ncol < V[t].cols(); ++ncol) {
        std::fill(eigen_vec.begin(), eigen_vec.end(), Complex(.0, .0));
        infile.read((char *)&(eigen_vec[0]), 2 * V[t].rows() * sizeof(double));
        if (!infile) {
          throw std::runtime_error("Problem while reading Eigenvectors!");
        }
        for (ssize_t nrow = 0; nrow < V[t].rows(); ++nrow) {
          (V[t])(nrow, ncol) = eigen_vec[nrow];
        }
      }
    } else {
      throw std::runtime_error("Eigenvector file does not exist!");
    }
    infile.close();

    // small test of trace and sum over the eigen vector matrix!
    if (verbose) {
      std::cout << "trace of V^d*V"
                << ":\t" << (V[t].adjoint() * V[t]).trace() << std::endl;
      std::cout << "sum over all entries of V^d*V"
                << ":\t" << (V[t].adjoint() * V[t]).sum() << std::endl;
    }
  } else {
    std::cout << "Randomizing eigenvectors ts: " << t << std::endl;
    V[t].setRandom();
  }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void EigenVector::read_eigen_vector(const std::string &filename, const ssize_t verbose) {
  for (int t = 0; t < ssize(V); t++) {
    std::string path = (boost::format("%s%03d") % filename % t).str();
    read_eigen_vector(path, t, verbose);
  }
}

void EigenVector::set_V(Eigen::MatrixXcd &v, const ssize_t t) {
  V[t] = v;
}
