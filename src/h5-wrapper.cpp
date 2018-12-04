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

template <>
H5::CompType make_comp_type<ComplexProduct>() {
  H5::CompType cmplxcmplx_w(4 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplxcmplx_w.insertMember("re_full", 0 * sizeof(double), type);
  cmplxcmplx_w.insertMember("im_full", 1 * sizeof(double), type);
  cmplxcmplx_w.insertMember("re_constrained", 2 * sizeof(double), type);
  cmplxcmplx_w.insertMember("im_constrained", 3 * sizeof(double), type);

  return cmplxcmplx_w;
}

#if 0
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
#endif

template <>
void write_homogenious(H5::DataSet &data_set,
                       std::vector<ComplexProduct> const &payload) {
  data_set.write(payload.data(), make_comp_type<ComplexProduct>());
}
