#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;

class LinearInterpolation{
  MatrixXd V;
  public:
    LinearInterpolation(const MatrixXd &V0){
      V = V0;
    }

    // evaluate function at time t
    float eval_function(float t){
      //we make an assumption that x are ordered
      float k = 0.0f;
      float b = 0.0f;
      for (int i = 0; i < V.rows()-1; i++){
        if ((t >= V(i,0)) && (t <= V(i+1,0))){
          k = (V(i,1) - V(i+1,1))/(V(i,0) - V(i+1,0));
          b = V(i,1) - V(i,0)*(V(i,1) - V(i+1,1))/(V(i,0) - V(i+1,0));
        }
      }
      return k*t + b;
      }
  };
