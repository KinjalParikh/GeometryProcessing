#include "poisson_surface_reconstruction.h"
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <igl/copyleft/marching_cubes.h>
#include <igl/cat.h>
#include <algorithm>
#include <igl/sum.h>

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, CONGRATULATIONS! You get this for free!
  ////////////////////////////////////////////////////////////////////////////
  // number of input points
  const int n = P.rows();
  // Grid dimensions
  int nx, ny, nz;
  // Maximum extent (side length of bounding box) of points
  double max_extent =
    (P.colwise().maxCoeff()-P.colwise().minCoeff()).maxCoeff();
  // padding: number of cells beyond bounding box of input points
  const double pad = 8;
  // choose grid spacing (h) so that shortest side gets 30+2*pad samples
  double h  = max_extent/double(30+2*pad);
  // Place bottom-left-front corner of grid at minimum of points minus padding
  Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
  // Grid dimensions should be at least 3 
  nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
  ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
  nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);
  // Compute positions of grid nodes
  Eigen::MatrixXd x(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++) 
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  //Find weights for staggered grid in each direction
  //std::cout << "Computing Wx \n";
  Eigen::SparseMatrix<double> Wx(n, ((nx-1) * ny * nz));
  fd_interpolate(nx-1, ny, nz, h, corner+Eigen::RowVector3d(h/2, 0, 0), P, Wx);

  //std::cout << "Computing Wy \n";
  Eigen::SparseMatrix<double> Wy(n, (nx * (ny-1) * nz));
  fd_interpolate(nx, ny-1, nz, h, corner + Eigen::RowVector3d(0, h / 2, 0), P, Wy);

  //std::cout << "Computing Wz \n";
  Eigen::SparseMatrix<double> Wz(n, (nx * ny * (nz-1)));
  fd_interpolate(nx, ny, nz-1, h, corner + Eigen::RowVector3d(0, 0, h / 2), P, Wz);

  //V  = W_transpose* N for each direction
  //std::cout << "Computing Vx \n";
  Eigen::MatrixXd Vx = Eigen::SparseMatrix<double>(Wx.transpose())* N.col(0);
  //std::cout << "Computing Vy \n";
  Eigen::MatrixXd Vy = Eigen::SparseMatrix<double>(Wy.transpose()) * N.col(1);
  //std::cout << "Computing Vz \n";
  Eigen::MatrixXd Vz = Eigen::SparseMatrix<double>(Wz.transpose()) * N.col(2);
  Eigen::MatrixXd Vcat(((nx - 1) * ny * nz) + (nx * (ny - 1) * nz) + (nx * ny * (nz - 1)), 1);
  igl::cat(1, igl::cat(1, Vx, Vy), Vz, Vcat);

  //Calculate G
  //std::cout << "Computing G \n";
  Eigen::SparseMatrix<double> G(((nx - 1) * ny * nz) + (nx * (ny - 1) * nz) + (nx * ny * (nz - 1)), nx * ny * nz);
  fd_grad(nx, ny, nz, h, G);


  // Gt G g  = Gt v 
  //std::cout << "Computing g \n";
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solverA;
  Eigen::SparseMatrix<double> A = G.transpose() * G;
  Eigen::MatrixXd b = G.transpose() * Vcat;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  solver.compute(A);
  g = solver.solve(b);
  /* ... update b ... */
  g = solver.solve(b); // solve again


  //std::cout << "Computing sigma \n";
  Eigen::SparseMatrix<double> W(n, nx * ny * nz);
  fd_interpolate(nx, ny, nz, h, corner, P, W); 
  double sigma = (Eigen::MatrixXd::Constant(1, n, 1) *(W*g)).sum()/n; 
  //std::cout << "sigma: " <<sigma << std::endl;
  g = g - Eigen::VectorXd::Constant(nx * ny * nz, sigma);
  
  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
