#include "EigenVector.hpp"
#include "GaugeField.hpp"
#include "typedefs.hpp"

#include <boost/format.hpp>

#include <string>

int main() {
  std::string path_ev_in("/hiskp4/bartek/eigensystems/quenched/wilson_b5.85_L12T24/data");
  std::string path_ev_out("/hiskp4/helmes/contractions/test/derivative_operators/covariant_shift");
  std::string path_gauge_in = std::string("/hiskp4/gauges/quenched/wilson_b5.85_L12T24");
  const ssize_t dim_row = 12*12*12*3;
  EigenVector V_t(24,dim_row,48);
  EigenVector UV(24,dim_row,48);
  EigenVector UVshift(24,dim_row,48);
  EigenVector derivative(24,dim_row,48);
  auto const path_in = (boost::format("%s/eigenvectors.0400.") %
                          path_ev_in).str();
  V_t.read_eigen_vector(path_in,0);
  // Read in timeslice of gauge
  GaugeField gauge = GaugeField(24,12,12,12, path_gauge_in,0,23,4);
  gauge.read_gauge_field(400,0,23);
  std::cout << V_t[1](0,0) << std::endl;
  std::cout << V_t[1].col(0).segment(3*(3*12*12+4*12+5),3) <<std::endl;

  //Test the displacement phase
  //DisplacementDirection displacements {std::pair<char,char> ('>','x'),
  //                             std::pair<char,char> ('>','y'),
  //                             std::pair<char,char> ('<','z')};
  DisplacementDirection displacements {std::pair<char,char> ('>','x'),
                               std::pair<char,char> ('>','y'),
                               std::pair<char,char> ('<','x'),
                               std::pair<char,char> ('<','y')
                               };
  for (auto& d : displacements) std::cout<<d.first<<' '<<d.second<<std::endl;
  Eigen::Vector3f displacement_phase_vector = gauge.summed_displacement(displacements);
  std::cout<< displacement_phase_vector<< std::endl;

  //Displace something
  for (ssize_t t = 0; t < 24; ++t){
    Eigen::MatrixXcd df_eigvecs = gauge.forward_uv(V_t[t],t,'x',1);
    Eigen::MatrixXcd db_eigvecs = gauge.backward_uv(V_t[t],t,'x',1);
//    Eigen::MatrixXcd disp_eigvecs = gauge.displace_eigenvectors(V_t[t],t,displacements,1);
 //   std::cout<<(V_t[t].adjoint()*disp_eigvecs).trace()<<std::endl;
    std::cout<<(V_t[t].adjoint()*V_t[t]).trace()<<std::endl;
    //// Calculate Umu_V
    //Eigen::MatrixXcd umu_v = gauge.Umu_times_V(V_t[t],t,0,1);
    //// Calculate Umu_V_shifted
    //Eigen::MatrixXcd umu_vshift = gauge.Umu_times_shiftedV(V_t[t],t,0,1);
    //// Calculate complete derivative
    //Eigen::MatrixXcd symmetric_derivative = gauge.symmetric_derivative(V_t[t],
    //                                                                   t,0);
    //std::cout << "Timeslice: " << t 
    //          << " tr(V^dagDV) = " << symmetric_derivative.trace() << std::endl;

    //// Write them to disk
    //UV.set_V(umu_v,t);
    //UVshift.set_V(umu_vshift,t);
    //derivative.set_V(symmetric_derivative,t);
    //auto const path_uxvx= (boost::format("%s/%s/eigenvectors.0400.%03d") 
    //                       % path_ev_out % "/UxVx" % t).str();
    //UV.write_eigen_vector(path_uxvx, t, 0);

    //auto const path_uxvxp1= (boost::format("%s/%s/eigenvectors.0400.%03d") 
    //                       % path_ev_out % "/UxVxp1" % t).str();
    //UVshift.write_eigen_vector(path_uxvxp1, t, 0);

    //auto const path_derivative= (boost::format("%s/%s/eigenvectors.0400.%03d") 
    //                       % path_ev_out % "/derivative" % t).str();
    //derivative.write_eigen_vector(path_derivative, t, 0);
  }
}
