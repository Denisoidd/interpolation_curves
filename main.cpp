#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include "lagrange.cpp"
#include "linear.cpp"
#include "cubic.cpp"
#include "hermite.cpp"

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V1; // matrix storing vertex coordinates of the input curve
MatrixXi F1;

/* Draw a curve given a list of points. Successive points in the list are connected by a segment*/
void draw_curve(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V){
  viewer.append_mesh();
  for (unsigned i = 0; i < V.rows()-1; ++i)
    viewer.data().add_edges(
        V.row(i),
        V.row(i+1),
        Eigen::RowVector3d(1, 0, 0));
}

/* Draw tangents on curve */
void draw_tangents(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V, const MatrixXd &dX){
  viewer.append_mesh();
  std::cout << "dX" << '\n';
  std::cout << dX << '\n';
  for (unsigned i = 0; i < V.rows()-2; ++i){
    float dx = dX(i*2+1,0);
    float dy = dX(i*2+1,1);
    float t = 0.3;
    std::cout << "dx and dy " << i+1 << '\n';
    std::cout << dx << " " << dy << '\n';
    viewer.data().add_edges(
        dX.row(i*2),
        Eigen::RowVector3d(dX(i*2,0) + t*dx,dX(i*2,1) + t*dy,0),
        Eigen::RowVector3d(1, 5, 0));
  }
}

/*draw points from the list of points V*/
void draw_points(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V){
  viewer.append_mesh();
  viewer.data(0).add_points(V, Eigen::RowVector3d(1, 0, 0));
}

/*interpolate the leftmost and rightmost X position*/
void build_linspace(MatrixXd &linspace,const MatrixXd &V){
  for (size_t i = 0; i < linspace.rows(); i++) {
    linspace(i, 0) = V.col(0).minCoeff() + ((V.col(0).maxCoeff() - V.col(0).minCoeff())/(linspace.rows() - 1))*i;
  }
}

// int main(int argc, char *argv[])
// {
//   igl::readPLY("../data/curve0.ply", V1, F1);
//   //  print the number of mesh elements
//   std::cout << "Points: " << V1.rows() << std::endl;
//
//   // choose interpolation method
//   //LinearInterpolation interp(V1);
//   //LagrangeInterpolation interp(V1);
//   CubicInterpolation interp(V1);
//
//   int resolution = 500; // resolution of displayed curve
//   MatrixXd linspace = MatrixXd::Zero(resolution, 3);
//   build_linspace(linspace, V1); // initialize the X axis of the interpolation
//
//   for (size_t i = 0; i < resolution; i++) {
//     linspace(i, 1) = interp.eval_function(linspace(i, 0));
//   }
//
//   igl::opengl::glfw::Viewer viewer; // create the 3d viewer
//   draw_points(viewer, V1); // draw the bounding box (red edges and vertices)
//   draw_curve(viewer, linspace);
//
//   // uncomment if we want to draw tangents
//   //    MatrixXd dX(V1.rows(), 2);
//   //   for (size_t i = 1; i < V1.rows()-1; i++) {
//   //     interp.eval_tangent(i, dX, V1(i, 0));
//   //   }
//   // draw_tangents(viewer, V1, dX);
//   viewer.launch(); // run the editor
// }


// use this main for exercice 3
int main(int argc, char *argv[])
{
  igl::readPLY("../data/curve0.ply", V1, F1);
  std::cout << "Points: " << V1.rows() << std::endl;
  HermiteInterpolation interp(V1);
  int resolution = 400;
  MatrixXd linspace = MatrixXd::Zero(resolution, 3);
  MatrixXd linspace_for_tangent = MatrixXd::Zero((V1.rows()-2)*2,3);
  interp.eval_function(linspace);
  interp.eval_tangent(linspace_for_tangent);
  // eval function updates the entirety of linspace
  // std::cout << "linspace" << '\n';
  // std::cout << linspace << '\n';
  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  draw_points(viewer, V1); // draw the bounding box (red edges and vertices)
  draw_curve(viewer, linspace);
  draw_tangents(viewer, V1, linspace_for_tangent);
  viewer.launch(); // run the editor
}
