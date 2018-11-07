#include "h5-wrapper.hpp"

#include "DilutedFactor.hpp"

#include <boost/filesystem.hpp>

void create_folder(std::string const &path) {
  if (access(path.c_str(), 0) != 0) {
    std::cout << "\tdirectory " << path << " does not exist and will be created";
    boost::filesystem::path dir(path.c_str());
    if (boost::filesystem::create_directories(dir))
      std::cout << "\tSuccess" << std::endl;
    else
      std::cout << "\tFailure" << std::endl;
  }
}

struct Complex1Times {
  int t;
  double re;
  double im;
};

template <>
H5::CompType make_comp_type<Complex>() {
  H5::CompType cmplx_w(2 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplx_w.insertMember("re", HOFFSET(complex_t, re), type);
  cmplx_w.insertMember("im", HOFFSET(complex_t, im), type);

  return cmplx_w;
}

template <>
H5::CompType make_comp_type<ComplexProduct>() {
  H5::CompType cmplxcmplx_w(4 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplxcmplx_w.insertMember("rere", HOFFSET(ComplexProduct, rere), type);
  cmplxcmplx_w.insertMember("reim", HOFFSET(ComplexProduct, reim), type);
  cmplxcmplx_w.insertMember("imre", HOFFSET(ComplexProduct, imre), type);
  cmplxcmplx_w.insertMember("imim", HOFFSET(ComplexProduct, imim), type);

  return cmplxcmplx_w;
}

template <>
H5::CompType make_comp_type<Complex1Times>() {
  H5::CompType comp_type(sizeof(Complex1Times));
  auto const type_int = H5::PredType::NATIVE_INT;
  auto const type = H5::PredType::NATIVE_DOUBLE;
  comp_type.insertMember("t", HOFFSET(Complex1Times, t), type_int);
  comp_type.insertMember("re", HOFFSET(Complex1Times, re), type);
  comp_type.insertMember("im", HOFFSET(Complex1Times, im), type);

  return comp_type;
}

template <>
void write_homogenious(H5::DataSet &data_set, std::vector<Complex> const &payload) {
  data_set.write(payload.data(), make_comp_type<Complex>());
}

template <>
void write_homogenious(H5::DataSet &data_set,
                       std::vector<ComplexProduct> const &payload) {
  data_set.write(payload.data(), make_comp_type<ComplexProduct>());
}

/** Creates compound datatype to write complex numbers from complex_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for complex numbers
 */
template <>
H5::CompType comp_type_factory<Complex>() {
  H5::CompType cmplx_w(2 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplx_w.insertMember("re", HOFFSET(complex_t, re), type);
  cmplx_w.insertMember("im", HOFFSET(complex_t, im), type);

  return cmplx_w;
}

/** Creates compound datatype to write complex numbers from compcomp_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for structs of four doubles
 */
template <>
H5::CompType comp_type_factory<ComplexProduct>() {
  H5::CompType cmplxcmplx_w(4 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplxcmplx_w.insertMember("rere", HOFFSET(ComplexProduct, rere), type);
  cmplxcmplx_w.insertMember("reim", HOFFSET(ComplexProduct, reim), type);
  cmplxcmplx_w.insertMember("imre", HOFFSET(ComplexProduct, imre), type);
  cmplxcmplx_w.insertMember("imim", HOFFSET(ComplexProduct, imim), type);

  return cmplxcmplx_w;
}
