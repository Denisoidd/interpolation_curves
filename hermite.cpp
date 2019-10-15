#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;
class HermiteInterpolation{

  MatrixXd solutionx;// solution for x
  MatrixXd solutiony;// solution for y
  MatrixXd steps;// points to interpolate
  MatrixXd solvey;// left size of linear system
  MatrixXd solvex;// left size of linear system
  MatrixXd slopes;
  MatrixXd W;

void solve_x(){
  int num_block = 0;
  for (int i = 1; i < steps.rows(); i++){
    MatrixXd part_solv_x = MatrixXd::Zero(4,1);

    solvex(0,0) = steps(i-1,0);
    solvex(1,0) = steps(i,0);
    solvex(2,0) = slopes(i-1,0);
    solvex(3,0) = slopes(i,0);

    part_solv_x = W.colPivHouseholderQr().solve(solvex);

    solutionx(num_block*4,0) = part_solv_x(0,0);
    solutionx(num_block*4+1,0) = part_solv_x(1,0);
    solutionx(num_block*4+2,0) = part_solv_x(2,0);
    solutionx(num_block*4+3,0) = part_solv_x(3,0);

    num_block++;
  }
}

void solve_y(){
  int num_block = 0;
  for (int i = 1; i < steps.rows(); i++){
    MatrixXd part_solv_y = MatrixXd::Zero(4,1);

    solvey(0,0) = steps(i-1,1);
    solvey(1,0) = steps(i,1);
    solvey(2,0) = slopes(i-1,1);
    solvey(3,0) = slopes(i,1);

    part_solv_y = W.colPivHouseholderQr().solve(solvey);

    solutiony(num_block*4,0) = part_solv_y(0,0);
    solutiony(num_block*4+1,0) = part_solv_y(1,0);
    solutiony(num_block*4+2,0) = part_solv_y(2,0);
    solutiony(num_block*4+3,0) = part_solv_y(3,0);

    num_block++;
  }
}

public:
  HermiteInterpolation(const MatrixXd &V1){
    W = MatrixXd::Zero(4,4);

    W(0,0) = 1;
    W(1,0) = 1;
    W(1,1) = 1;
    W(1,2) = 1;
    W(1,3) = 1;
    W(2,1) = 1;
    W(3,1) = 1;
    W(3,2) = 2;
    W(3,3) = 3;

    steps = V1;

    slopes = MatrixXd::Zero(5,2);
    slopes(0,0) = 1; slopes(0,1) = 1;
    slopes(1,0) = 1; slopes(1,1) = 1;
    slopes(2,0) = 1.0; slopes(2,1) = 1;
    slopes(3,0) = -2.5; slopes(3,1) = 1;
    slopes(4,0) = 1.0; slopes(4,1) = 1;
    solvex = MatrixXd::Zero(4,1);
    solvey = MatrixXd::Zero(4,1);
    solutionx = MatrixXd::Zero((V1.rows()-1)*4,1);
    solutiony = MatrixXd::Zero((V1.rows()-1)*4,1);
    solve_x();
    solve_y();
  }

  // complete linspace with corrdinates
  void eval_function(MatrixXd &linspace){

    float d_t_x = (1.0f/linspace.rows())*(steps.rows()-1);
    float d_t_y = (1.0f/linspace.rows())*(steps.rows()-1);

    float t_x = 0;
    float t_y = 0;

    for (int k = 0; k < steps.rows()-1; k++){
      float t_x = 0;
      float t_y = 0;
      for (int i = (k)*(linspace.rows()/(steps.rows()-1)); i < (k+1)*(linspace.rows()/(steps.rows()-1)); i++){
        linspace(i,0) = solutionx(k*4,0) + solutionx(k*4+1,0)*t_x + solutionx(k*4+2,0)*std::pow(t_x, 2) + solutionx(k*4+3,0)*std::pow(t_x, 3);
        linspace(i,1) = solutiony(k*4,0) + solutiony(k*4+1,0)*t_y + solutiony(k*4+2,0)*std::pow(t_y, 2) + solutiony(k*4+3,0)*std::pow(t_y, 3);
        t_x = t_x + d_t_x;
        t_y = t_y + d_t_y;
      }
    }
    t_x = 0;
    t_y = 0;
      for (int i = (steps.rows() - 2) * (linspace.rows()/(steps.rows()-1)); i < linspace.rows(); i++){
        linspace(i,0) = solutionx((steps.rows() - 2)*4,0) + solutionx((steps.rows() - 2)*4+1,0)*t_x + solutionx((steps.rows() - 2)*4+2,0)*std::pow(t_x, 2) + solutionx((steps.rows() - 2)*4+3,0)*std::pow(t_x, 3);
        linspace(i,1) = solutiony((steps.rows() - 2)*4,0) + solutiony((steps.rows() - 2)*4+1,0)*t_y + solutiony((steps.rows() - 2)*4+2,0)*std::pow(t_y, 2) + solutiony((steps.rows() - 2)*4+3,0)*std::pow(t_y, 3);
        t_x = t_x + d_t_x;
        t_y = t_y + d_t_y;
      }
  }

  void eval_tangent(MatrixXd &linspace){
    for (int i=0; i < steps.rows()-2; i++){

      linspace.row(i*2) = steps.row(i+1);

      linspace(i*2+1,0) = solutionx(i*4+1,0) + 2*solutionx(i*4+2,0) + 3*solutionx(i*4+3,0);
      linspace(i*2+1,1) = solutiony(i*4+1,0) + 2*solutiony(i*4+2,0) + 3*solutiony(i*4+3,0);

    }
    std::cout << "Linspace" << '\n';
    std::cout << linspace << '\n';
  }
};
