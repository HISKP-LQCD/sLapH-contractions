#include "Perambulator.hpp"

#include "StopWatch.hpp"

/**
 *  @param entity       The entry where this peram will be stored
 *  @param Lt           Total number of timeslices - for each peram the same
 *  @param nb_eigen_vec Total number of eigen vecs - for each peram the same
 *  @param quark        Contains information about dilution scheme and size
 *  @param filename     Just the file name
 */
void Perambulator::read_perambulator(const ssize_t entity,
                                     const ssize_t Lt,
                                     const ssize_t nb_eigen_vec,
                                     const quark &quark,
                                     const std::string &filename,
                                     const bool mock) {
  StopWatch swatch("Perambulator I/O");
  swatch.start();
  if (!mock) {
    FILE *fp = NULL;
    std::cout << "\tReading perambulator from file:\n\t\t" << filename << "\n";

    // reading the data into temporary array
    std::vector<Complex> perambulator_read(peram[entity].size());
    if ((fp = fopen(filename.c_str(), "rb")) == NULL) {
      std::cout << "failed to open file to read perambulator: " << filename << "\n"
                << std::endl;
      exit(1);
    }
    int check_read =
        fread(&(perambulator_read[0]), sizeof(Complex), peram[entity].size(), fp);
    fclose(fp);
    // check if all data were read in
    if (check_read != peram[entity].size()) {
      std::cout << "\n\nFailed to read perambulator. Expected size "
                << peram[entity].size() << " and got " << check_read << " B."
                << std::endl;
      exit(1);
    }

    // setting peram to zero
    peram[entity].setZero();

    // re-sorting and copy into matrix structure
    // TODO: At this point it is very easy to included different dilution schemes.
    //       However, due to simplicity this will be postponed!
    const ssize_t nb_dil_T = quark.number_of_dilution_T;
    const ssize_t nb_dil_E = quark.number_of_dilution_E;
    const ssize_t nb_dil_D = quark.number_of_dilution_D;
    ssize_t col_i, row_i;
    const int nb_inversions = Lt * nb_dil_E * nb_dil_D / nb_dil_T;
    for (ssize_t t1 = 0; t1 < Lt; ++t1)
      for (ssize_t ev1 = 0; ev1 < nb_eigen_vec; ++ev1)
        for (ssize_t dirac1 = 0; dirac1 < 4; ++dirac1)
          for (ssize_t t2 = 0; t2 < (Lt / nb_dil_T); ++t2)
            for (ssize_t ev2 = 0; ev2 < nb_dil_E; ++ev2)
              for (ssize_t dirac2 = 0; dirac2 < nb_dil_D; ++dirac2) {
                row_i = 4 * nb_eigen_vec * t1 + 4 * ev1 + dirac1;
                col_i = nb_dil_D * nb_dil_E * t2 + nb_dil_D * ev2 + dirac2;
                peram[entity](4 * nb_eigen_vec * t1 + nb_eigen_vec * dirac1 + ev1,
                              nb_dil_E * nb_dil_D * t2 + nb_dil_E * dirac2 + ev2) =
                    perambulator_read[row_i * nb_inversions + col_i];
              }
  } else {
    std::cout << "\tRandomizing perambulator" << std::endl;
    peram[entity].setRandom();
  }
  swatch.stop();
  swatch.print();
}

/**
 *  @param Lt            Total number of timeslices - for each peram the same
 *  @param nb_eigen_vec  Total number of eigen vecs - for each peram the same
 *  @param quark         Contains information about dilution scheme and size
 *  @param filename_list Vector which contains all file names
 *
 *  Loops over \code nb_entities = quark.size() * nb_rnd_vec \endcode. For each
 *  internally read_perambulator() for a single file is called
 */
void Perambulator::read_perambulators_from_separate_files(
    const ssize_t Lt,
    const ssize_t nb_eigen_vec,
    const std::vector<quark> &quark,
    const std::vector<std::string> &filename_list) {
  if (filename_list.size() != peram.size())
    std::cout << "Problem when reading perambulators: The number of "
              << "perambulators read is not the same as the expected one!" << std::endl;
  ssize_t j = 0;  // TODO: Not beautiful but practical - Think about change
  for (ssize_t i = 0; i < ssize(quark); i++) {
    for (ssize_t r = 0; r < quark[i].number_of_rnd_vec; r++) {
      read_perambulator(j, Lt, nb_eigen_vec, quark[i], filename_list[j]);
      j++;
    }
  }
}
