#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;

class LagrangeInterpolation{
  MatrixXd W;
  MatrixXd a;

  MatrixXd Vandermonde(const MatrixXd &V, MatrixXd &W){
    double cur_elem = 0.0d;
    //VectorXd iden_vec = VectorXd::Ones(pol_deg + 1);
    //std::cout << iden_vec << std::endl;
    for (int i=0; i < W.rows(); i++){
      cur_elem = V(i,0);
      for (int j=0; j < W.cols(); j++){
        W(i,j) = std::pow(cur_elem, j);
      }
    }
    return W;
  }

public:
  LagrangeInterpolation(const MatrixXd &V1){
    W = Eigen::MatrixXd::Ones(V1.rows(),V1.rows());
    W = Vandermonde(V1, W);
    a = W.colPivHouseholderQr().solve(V1.col(1));
  }

  // evaluate interpolated function at time t
  float eval_function(float t){
    float res = a(0,0);
    for (int i=1; i < a.rows(); i++){
      res = res + std::pow(t,i) * a(i, 0);
    }
    return res;
  }
};
