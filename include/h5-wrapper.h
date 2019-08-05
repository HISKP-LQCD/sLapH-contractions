#pragma once

#include "typedefs.h"

#include "H5Cpp.h"

#include <boost/filesystem.hpp>

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include <omp.h>

#define PRINT(x) std::cout << #x << ": " << (x) << std::endl;

void create_folder(std::string const &path);

class HDF5Handle {
 public:
  HDF5Handle(std::string const &output_path,
             std::string const &diagram,
             std::string const &output_filename) {
    H5::Exception::dontPrint();

    create_folder(output_path);
    H5std_string const file_name((output_path + "/" + diagram + output_filename).c_str());
    file_ = H5::H5File(file_name, H5F_ACC_TRUNC);
  }

  H5::DataSet create_data_set(std::string const &dataset_name,
                              H5::CompType const &comp_type,
                              H5::DataSpace const &data_space) {
    return file_.createDataSet(dataset_name.c_str(), comp_type, data_space);
  }

  H5::Group create_group(std::string const &name) {
    return file_.createGroup(name.c_str());
  }

 private:
  H5::H5File file_;
};

template <typename Payload>
void write_homogenious(H5::DataSet &handle, Payload const &payload);

template <typename Payload>
void write_heterogenious(H5::Group &handle, Payload const &payload);

template <typename T>
H5::CompType make_comp_type();

template <size_t rvecs>
void write_heterogenious(
    H5::Group &group,
    boost::detail::multi_array::sub_array<std::vector<DilutedTrace<rvecs>>, 1ul> const
        &payload) {
  std::map<std::string, std::vector<Complex>> data;

  for (int t = 0; t < payload.size(); ++t) {
    PRINT(payload[t].size());
    auto const &subs = sub_accumulate(payload[t]);
    PRINT(subs.size());

    for (auto const &elem : subs) {
      std::ostringstream oss;
      bool first = true;
      for (auto const &id : elem.first) {
        if (first) {
          first = false;
        } else {
          oss << "_";
        }

        oss << "rnds";
        // The `static_cast<int>` is needed because otherwise the `uint8_t` is interpreted
        // as `signed char` and then this will define the string representation. We want
        // it interpreted as an integer.
        oss << static_cast<int>(id);
      }
      std::string const key = oss.str();

      std::cout << key << " " << elem.second << std::endl;

      data[key].push_back(elem.second);
    }
  }

  auto const &comp_type = make_comp_type<Complex>();

  for (auto const &pair : data) {
    hsize_t dim(pair.second.size());
    H5::DataSpace dspace(1, &dim);
    auto data_set = group.createDataSet(pair.first.c_str(), comp_type, dspace);
    data_set.write(pair.second.data(), comp_type);
  }
}

template <typename Numeric>
H5::CompType comp_type_factory();

/*! Class to write correlations function to files in hdf5 format
 *
 *  @Warning  Dependency inversion principle is violated: Class depends on the
 *            concrete implementation of hdf5
 */
class WriteHDF5Correlator {
 public:
  WriteHDF5Correlator(const std::string output_path,
                      const std::string diagram,
                      const std::string output_filename,
                      const H5::CompType &_comp_type)
      : comp_type(_comp_type) {
    create_folder_for_hdf5_file(output_path.c_str());

    std::string const file_name(output_path + "/" + diagram + output_filename);
    open_or_create_hdf5_file(file_name);
  }

  /*! Writes data to file
   *
   *  @Param corr       The data to write
   *  @Param corr_info  Contains the hdf5_dataset_name
   *
   *  @todo   It is sufficient to pass the hdf5_dataset_name
   *  @todo   corr_datatype would be better as class template typename
   *
   *  @remark The type corr_datatype is always either complex_t or
   *          compcomp_t. The function body is identical for both types
   *          as everything is specified by corr_info. Thus the template
   *          overload
   */
  template <typename corr_datatype>
  void write(const std::vector<corr_datatype> &corr, const DiagramIndex &corr_info) {
    // Exceptions are automatically printed, we do not need this feature.
    H5::Exception::dontPrint();

    // Create a data set object.
    H5::Group group;
    H5std_string dataset_name((corr_info.hdf5_dataset_name).c_str());

    hsize_t dim(corr.size());
    H5::DataSpace dspace(1, &dim);

    // We try to open the data set in the file. If it exists, we do not need to do
    // anything because we do not overwrite existing data.
    try {
      file.openDataSet(dataset_name);
      std::cout << "Not writing " << corr_info.hdf5_dataset_name << " because it exists."
                << std::endl;
      return;
    } catch (H5::Exception &) {
    }

    // Actual write.
    try {
      auto dset = file.createDataSet(dataset_name, comp_type, dspace);
      dset.write(&corr[0], comp_type);
    } catch (H5::Exception &e) {
      e.printError();
    }
  }

 private:
  /*! Checks whether output path exists and if not creates it
   *
   *  @param[in] path Path where hdf5 file shall be written
   */
  void create_folder_for_hdf5_file(const char *path) {
    if (access(path, 0) != 0) {
      std::cout << "\tdirectory " << path << " does not exist and will be created";
      boost::filesystem::path dir(path);
      if (boost::filesystem::create_directories(dir))
        std::cout << "\tSuccess" << std::endl;
      else
        std::cout << "\tFailure" << std::endl;
    }
  }

  /*! Opens correlator file, preserving existing data.
   *
   *  @param[in] name String containing path+filename of the desired file
   */
  void open_or_create_hdf5_file(std::string const &name) {
    const H5std_string file_name(name.c_str());

    if (boost::filesystem::exists(name)) {
      file = H5::H5File(name, H5F_ACC_RDWR);
    } else{
      file = H5::H5File(name, H5F_ACC_TRUNC);
    }
  }

  /*! The hdf5 file pointer */
  H5::H5File file;
  /*! The hdf5 compound datatype.
   *
   *  @see  H5::CompType comp_type_factory_tr()
   *  @see  H5::CompType comp_type_factory<compcomp_t>()
   */
  H5::CompType comp_type;

};  // end of class WriteHDF5Correlator
