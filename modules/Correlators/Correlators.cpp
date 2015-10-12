#include "Correlators.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
LapH::Correlators::Correlators (
                        const size_t Lt, const size_t dilT, const size_t dilE, 
                        const size_t nev, const CorrelatorLookup& corr_lookup) :
                                      Lt(Lt), dilT(dilT), dilE(dilE), nev(nev) {

  C1.resize(boost::extents[corr_lookup.C1.size()][Lt]);
  C20.resize(boost::extents[corr_lookup.C20.size()][Lt]);
  std::fill(C20.data(), C20.data()+C20.num_elements(), cmplx(.0,.0));
  C2c.resize(boost::extents[corr_lookup.C2c.size()][Lt]);
  std::fill(C2c.data(), C2c.data()+C2c.num_elements(), cmplx(.0,.0));

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C1(const Quarklines& quarklines,
                     const std::vector<CorrInfo>& corr_lookup,
                     const QuarklineLookup& quark_lookup,
                     const std::vector<RandomIndexCombinationsQ1>& ric_lookup) {

  std::cout << "\tcomputing C1:";
  clock_t time = clock();

  for(const auto& c_look : corr_lookup){
//    const auto& ric0 = ric_lookup[quark_lookup.Q1[c_look.lookup[0]].
//                                                     id_ric_lookup].rnd_vec_ids;
//    const auto& ric1 = ric_lookup[quark_lookup.Q1[c_look.lookup[1]].
//                                                     id_ric_lookup].rnd_vec_ids;
//    if(ric0.size() != ric1.size()){
//      std::cout << "rnd combinations are not the same in build_corr0" 
//                << std::endl;
//      exit(0);
//    }
//    for(size_t t1 = 0; t1 < Lt; t1++){
//    for(size_t t2 = 0; t2 < Lt; t2++){
//      corr0[c_look.id][t1][t2].resize(ric0.size());
//      for(auto& corr : corr0[c_look.id][t1][t2])
//        corr = cmplx(0.0,0.0);
//      for(const auto& rnd : ric0){
//        const auto id = &rnd - &ric0[0];
//        const auto it1 = std::find_if(ric1.begin(), ric1.end(),
//                                [&](std::pair<size_t, size_t> pair){
//                                  return (pair == 
//                                         std::make_pair(rnd.second, rnd.first));
//                                });
//        if(it1 == ric1.end()){
//          std::cout << "something wrong with random vectors in build_corr0" 
//                    << std::endl;
//          exit(0);
//        }
//        //for(size_t block = 0; block < 4; block++)
//          corr0[c_look.id][t1][t2][id] += 
//                      (quarklines.return_Q1(t1, t2/dilT, c_look.lookup[0], id) *
//                       quarklines.return_Q1(t2, t1/dilT, c_look.lookup[1], 
//                                                     it1-ric1.begin())).trace();
//      }
//    }}
  }
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corr0(const Quarklines& quarklines,
                     const std::vector<CorrInfo>& corr_lookup,
                     const QuarklineLookup& quark_lookup,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup) {

  std::cout << "\tcomputing corr0:";
  clock_t time = clock();

  corr0.resize(boost::extents[corr_lookup.size()][Lt][Lt]);
  for(const auto& c_look : corr_lookup){
    const auto& ric0 = ric_lookup[quark_lookup.Q1[c_look.lookup[0]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = ric_lookup[quark_lookup.Q1[c_look.lookup[1]].
                                                     id_ric_lookup].rnd_vec_ids;
    if(ric0.size() != ric1.size()){
      std::cout << "rnd combinations are not the same in build_corr0" 
                << std::endl;
      exit(0);
    }
    for(size_t t1 = 0; t1 < Lt; t1++){
    for(size_t t2 = 0; t2 < Lt; t2++){
      corr0[c_look.id][t1][t2].resize(ric0.size());
      for(auto& corr : corr0[c_look.id][t1][t2])
        corr = cmplx(0.0,0.0);
      for(const auto& rnd : ric0){
        const auto id = &rnd - &ric0[0];
        const auto it1 = std::find_if(ric1.begin(), ric1.end(),
                                [&](std::pair<size_t, size_t> pair){
                                  return (pair == 
                                         std::make_pair(rnd.second, rnd.first));
                                });
        if(it1 == ric1.end()){
          std::cout << "something wrong with random vectors in build_corr0" 
                    << std::endl;
          exit(0);
        }
        corr0[c_look.id][t1][t2][id] += 
                    (quarklines.return_Q1(t1, t2/dilT, c_look.lookup[0], id) *
                     quarklines.return_Q1(t2, t1/dilT, c_look.lookup[1], 
                                                     it1-ric1.begin())).trace();
      }
    }}
  }
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C20(const std::vector<CorrInfo>& corr_lookup) {


  for(const auto& c_look : corr_lookup){
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t1 - t2 - (int)Lt) % (int)Lt);
      for(const auto& corr : corr0[c_look.lookup[0]][t1][t2])
        C20[c_look.id][t] += corr;
    }}
    for(auto& corr : C20[c_look.id]){
      // normalisation
      corr /= Lt*corr0[c_look.id][0][0].size();
      std::cout << std::setprecision(5) << &corr - &C20[c_look.id][0] << "\t" 
                << corr << std::endl;
    }
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corrC(const Quarklines& quarklines, 
                     const OperatorsForMesons& meson_operator,
                     const OperatorLookup& operator_lookup,
                     const std::vector<CorrInfo>& corr_lookup,
                     const QuarklineLookup& quark_lookup) {

  std::cout << "\tcomputing corrC:";
  clock_t time = clock();

  corrC.resize(boost::extents[corr_lookup.size()][Lt][Lt]);
  for(const auto& c_look : corr_lookup){
    const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2V[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[//needed only for checking
               operator_lookup.rvdaggervr_lookuptable[c_look.lookup[1]].
                                                    id_ricQ_lookup].rnd_vec_ids;
    if(ric0.size() != ric1.size()){
      std::cout << "rnd combinations are not the same in build_corrC" 
                << std::endl;
      exit(0);
    }
    for(size_t t1 = 0; t1 < Lt; t1++){
    for(size_t t2 = 0; t2 < Lt; t2++){
      corrC[c_look.id][t1][t2].resize(ric0.size());
      for(auto& corr : corrC[c_look.id][t1][t2])
        corr = cmplx(0.0,0.0);
      for(const auto& rnd : ric0){
        const auto id = &rnd - &ric0[0];
        for(size_t block = 0; block < 4; block++)
          corrC[c_look.id][t1][t2][id] += 
               quarklines.return_gamma_val(c_look.gamma[0], block) *
               (quarklines.return_Q2V(t1, t2/dilT, c_look.lookup[0], id).
                                  block(block*dilE, block*dilE, dilE, dilE) *
               meson_operator.return_rvdaggervr(c_look.lookup[1], t2, id).block(
                       quarklines.return_gamma_row(c_look.gamma[0], block)*dilE,
                                               block*dilE, dilE, dilE)).trace();
      }
    }}
  }
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C2c(const std::vector<CorrInfo>& corr_lookup) {

  for(const auto& c_look : corr_lookup){
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t1 - t2 - (int)Lt) % (int)Lt);
      for(const auto& corr : corrC[c_look.lookup[0]][t1][t2]){
        C2c[c_look.id][t] += corr;
      }
    }}
    for(auto& corr : C2c[c_look.id]){
      // normalisation
      corr /= Lt*corrC[c_look.id][0][0].size();
      std::cout << std::setprecision(5) << &corr - &C2c[c_look.id][0] << "\t" 
                << corr << std::endl;
    }
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::contract (const Quarklines& quarklines, 
                     const OperatorsForMesons& meson_operator,
                     const OperatorLookup& operator_lookup,
                     const CorrelatorLookup& corr_lookup, 
                     const QuarklineLookup& quark_lookup) {

  // 1. Build all functions which need corrC and free it afterwards.
  build_corrC(quarklines, meson_operator, operator_lookup, corr_lookup.corrC, 
                                                                  quark_lookup);
  build_C2c(corr_lookup.C2c);
//  corrC.resize(boost::extents[0][0][0][0][0][0]);
  // 2. Build all functions which need corr0 and free it afterwards.
  build_corr0(quarklines, corr_lookup.corr0, quark_lookup, 
                                                  operator_lookup.ricQ2_lookup);
  build_C20(corr_lookup.C20);
//  corr0.resize(boost::extents[0][0][0][0][0][0]);
  // 3. Build all other correlation functions.
  build_C1(quarklines, corr_lookup.C1, quark_lookup, 
                                                  operator_lookup.ricQ1_lookup);
}





