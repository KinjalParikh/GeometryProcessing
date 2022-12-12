#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include "igl/cat.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  Eigen::SparseMatrix<double> Dx ((nx - 1) * ny * nz, nx * ny * nz);
	fd_partial_derivative(nx, ny, nz, h, 0, Dx);

	Eigen::SparseMatrix<double> Dy(nx * (ny-1) * nz, nx * ny * nz);
	fd_partial_derivative(nx, ny, nz, h, 1, Dy);

	Eigen::SparseMatrix<double> Dz(nx * ny * (nz - 1), nx * ny * nz);
	fd_partial_derivative(nx, ny, nz, h, 2, Dz);

	Eigen::SparseMatrix<double> G1((nx - 1) * ny * nz + nx * (ny - 1) * nz , nx * ny * nz);
	igl::cat(1, Dx, Dy, G1);
	igl::cat(1, G1, Dz, G);
}
