#include "global_data.hpp"

#include <boost/format.hpp>

namespace {

static std::vector<std::string> create_rnd_vector_file_names(
    const int config, const int nb_of_eigen_vec, const std::vector<quark> quarks) {
  std::vector<std::string> filename_list;  // the return vector of file names
  // running over all quarks
  for (const auto &q : quarks) {
    for (int rnd_vec_i = 0; rnd_vec_i < q.number_of_rnd_vec; ++rnd_vec_i) {
      // building paths and filenames for rnd vecs
      auto const path =
          (boost::format(
               "%s/cnfg%04d/rnd_vec_%02d/randomvector.rndvecnb%02d.%s.nbev%04d.%04d") %
           q.path % config % rnd_vec_i % rnd_vec_i % q.type.c_str() % nb_of_eigen_vec %
           config)
              .str();
      filename_list.push_back(path);
    }
  }

  return filename_list;
}

static std::vector<std::string> create_perambulator_file_names(
    const int config,
    const int Lt,
    const int Lx,
    const int Ly,
    const int Lz,
    const std::vector<quark> quarks) {
  std::vector<std::string> filename_list;  // the return vector of file names

  // running over all quarks
  for (const auto &q : quarks) {
    char dil_scheme_T = q.dilution_T.back();
    char dil_scheme_E = q.dilution_E.back();
    char dil_scheme_D = q.dilution_D.back();

    for (int rnd_vec_i = 0; rnd_vec_i < q.number_of_rnd_vec; ++rnd_vec_i) {
      // data path for qbig contractions
      auto const path =
          (boost::format("%s/cnfg%04d/rnd_vec_%02d/"
                         "perambulator.rndvecnb%02d.%s.Tso%c%04d.Vso%c%04d.Dso%c%1d.TsiF%"
                         "04d.SsiF%04d.DsiF4.CsiF3.smeared0.%05d") %
           q.path % config % rnd_vec_i % rnd_vec_i % q.type.c_str() % dil_scheme_T %
           (Lt / q.number_of_dilution_T) % dil_scheme_E % q.number_of_dilution_E %
           dil_scheme_D % q.number_of_dilution_D % Lt % (Lx * Ly * Lz) % config)
              .str();
      filename_list.push_back(path);
    }
  }
  return filename_list;
}

static std::string create_eigenvector_file_name(const ssize_t config,
                                                const std::string &path_eigenvectors,
                                                const std::string &name_eigenvectors) {
  return (boost::format("%s/%s.%04d.") % path_eigenvectors % name_eigenvectors %
          static_cast<int>(config))
      .str();
}

static std::string create_correlator_file_name(const ssize_t config) {
  return (boost::format("_cnfg%04d.h5") % static_cast<int>(config)).str();
}

}  // end anonymous namespace

/*!
 * @param config Configuration number for which paths and file names shall be
 *               constructed
 */
void build_IO_names(GlobalData &gd, const ssize_t config) {
  gd.rnd_vec_construct.filename_list =
      create_rnd_vector_file_names(config, gd.number_of_eigen_vec, gd.quarks);
  gd.peram_construct.filename_list =
      create_perambulator_file_names(config, gd.Lt, gd.Lx, gd.Ly, gd.Lz, gd.quarks);
  gd.filename_eigenvectors =
      create_eigenvector_file_name(config, gd.path_eigenvectors, gd.name_eigenvectors);

  gd.filename_ending_correlators = create_correlator_file_name(config);
}
