#include "fd_partial_derivative.h"
#include<iostream>

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  double val = 1 / h;
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;

	if (dir == 0) {
		tripletList.reserve((nx-1)*ny*nz);
		for (int i = 0; i < nx - 1; i++) {
			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < nz; k++) {
					int indexX = i + j * (nx-1) + k * ny * (nx-1);
					int indexY = i + j * nx + k * ny * nx;
					tripletList.push_back(T(indexX, indexY, -val));
					tripletList.push_back(T(indexX, indexY+1, val));
					if (indexX >= (nx - 1) * ny * nz) { 
						std::cout << "Out of range row index: " << i << " " << j << " " << k << std::endl; 
					}
					if (indexY+1>= nx * ny * nz) { 
						std::cout << "Out of range column index: " << i << " " << j << " " << k << std::endl; 
					}
					//D.insert(index, index) = -val;
					//D.insert(index, index + 1) = val;
				}
			}
		}
	}

	if (dir == 1) {
		tripletList.reserve((ny - 1) * nx * nz);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny - 1; j++) {
				for (int k = 0; k < nz; k++) {
					int indexX = i + j * nx + k * (ny - 1) * nx;
					int indexY = i + j * nx + k * ny * nx;
					tripletList.push_back(T(indexX, indexY, -val));
					tripletList.push_back(T(indexX, indexY + nx, val));
				}
			}
		}
	}

	if (dir == 2) {
		tripletList.reserve((nz - 1) * ny * nx);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < nz-1 ; k++) {
					int index = i + j * nx + k * ny * nx;
					tripletList.push_back(T(index, index, -val));
					tripletList.push_back(T(index, index + nx*ny, val));
				}
			}
		}
	}

	D.setFromTriplets(tripletList.begin(), tripletList.end());
}
