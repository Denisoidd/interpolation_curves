#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;

class CubicInterpolation{
  MatrixXd M;// coefficients
  MatrixXd a;// solution
  MatrixXd V;// points to interpolate
  MatrixXd y;// left size of linear system

  /*
  Initialize system constraints
  */
  void init_system(const MatrixXd &V, MatrixXd &M, MatrixXd &y){

    int r = 0;

    //std::cout << V << std::endl;
    for (int i=1; i < y.rows()-1; i++){
      int k = i % 4;
      if (k == 0){
        y(i,0) = 0;
        //std::cout << "1 / 4" << std::endl;
      }
      else if (k == 1){
        y(i,0) = V(r,1);
        //std::cout << "2 / 4" << std::endl;
        r = r + 1;
        //std::cout << "r" << std::endl;
        //std::cout << r << std::endl;
      }
      else if (k == 2){
        y(i,0) = V(r,1);
        //std::cout << "3 / 4" << std::endl;
        //std::cout << "r" << std::endl;
        //std::cout << r << std::endl;

      }
      else {
        y(i,0) = 0;
        //std::cout << "4 / 4" << std::endl;
      }
    }
    //std::cout << y << std::endl;
  }

public:
  CubicInterpolation(const MatrixXd &V1){
    M = MatrixXd::Zero(4*(V1.rows() - 1),4*(V1.rows() - 1));
    y = MatrixXd::Zero(4*(V1.rows() - 1),1);
    V = V1;
    init_system(V1, M, y );
    //first line
    M(0,1) = 1;
    M(0,2) = 2 * V(0,0);
    M(0,3) = 3 * V(0,0) * V(0,0);

    //all lines in the matrix
    int num_block = 0;
    for (int i = 1; i < V.rows() - 1; i++){

      M(num_block*4 + 1,num_block*4 + 0) = 1;
      M(num_block*4 + 1,num_block*4 + 1) = V(i-1,0);
      M(num_block*4 + 1,num_block*4 + 2) = std::pow(V(i-1,0),2);
      M(num_block*4 + 1,num_block*4 + 3) = std::pow(V(i-1,0),3);

      M(num_block*4 + 2,num_block*4 + 0) = 1;
      M(num_block*4 + 2,num_block*4 + 1) = V(i,0);
      M(num_block*4 + 2,num_block*4 + 2) = std::pow(V(i,0),2);
      M(num_block*4 + 2,num_block*4 + 3) = std::pow(V(i,0),3);

      M(num_block*4 + 3,num_block*4 + 1) = 1;
      M(num_block*4 + 3,num_block*4 + 2) = 2 * V(i,0);
      M(num_block*4 + 3,num_block*4 + 3) = 3 * std::pow(V(i,0),2);
      M(num_block*4 + 3,num_block*4 + 5) = -1;
      M(num_block*4 + 3,num_block*4 + 6) = -2 * V(i,0);
      M(num_block*4 + 3,num_block*4 + 7) = -3 * std::pow(V(i,0), 2);

      M(num_block*4 + 4,num_block*4 + 2) = 2;
      M(num_block*4 + 4,num_block*4 + 3) = 6 * V(i,0);
      M(num_block*4 + 4,num_block*4 + 6) = -2;
      M(num_block*4 + 4,num_block*4 + 7) = -6 * V(i,0);

      num_block++;
    }

    //the last line
    std::cout << "Num_block number" << '\n';
    std::cout << num_block << '\n';
    M(num_block*4+1,num_block*4) = 1;
    M(num_block*4+1,num_block*4+1) = V(V.rows() - 2, 0);
    M(num_block*4+1,num_block*4+2) = std::pow(V(V.rows() - 2, 0),2);
    M(num_block*4+1,num_block*4+3) = std::pow(V(V.rows() - 2, 0),3);

    M(num_block*4+2,num_block*4) = 1;
    M(num_block*4+2,num_block*4+1) = V(V.rows() - 1, 0);
    M(num_block*4+2,num_block*4+2) = std::pow(V(V.rows() - 1, 0),2);
    M(num_block*4+2,num_block*4+3) = std::pow(V(V.rows() - 1, 0),3);

    M(num_block*4+3,num_block*4+1) = 1;
    M(num_block*4+3,num_block*4+2) = 2 * V(V.rows() - 1, 0);
    M(num_block*4+3,num_block*4+3) = 3* std::pow(V(V.rows() - 1, 0),2);










    //M(num_block*4+3,num_block*4) = 1;
    // M(num_block*4+3,num_block*4+1) = 1;
    // M(num_block*4+3,num_block*4+2) = 2*V(V.rows() - 1, 0);
    // M(num_block*4+3,num_block*4+3) = 3*std::pow(V(V.rows() - 1, 0),2);
    std::cout << "Matrix M" << '\n';
    std::cout << M << '\n';
    std::cout << "vector y" << '\n';
    std::cout << y << '\n';
    std::cout << "Matrix V" << '\n';
    std::cout << V << '\n';

    a = M.colPivHouseholderQr().solve(y);
  }


  /*
  Evaluate tangent at step i
  */
  void eval_tangent(float i, MatrixXd &dX, float x){
    // complete here
  }

  /*
  Evaluate function at time t
  */
  float eval_function(float t){
    for (int i = 0; i < V.rows()-1; i++){
      if ((t >= V(i,0)) && (t <= V(i+1,0))){
        return a(i*4,0) + a(i*4 + 1,0)*t + a(i*4 + 2,0) * std::pow(t,2) + a(i*4 + 3,0) * std::pow(t,3);
      }
    }
  }
};
