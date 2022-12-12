#include "fd_interpolate.h"
#include <iostream>
#include <Eigen/Dense>

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  W.reserve(Eigen::VectorXd::Constant(nx*ny*nz, 4));
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(P.rows()*8);

	for (int p = 0; p < P.rows(); p++) {

		//if (p == 0) { std::cout << "P(" << p << "): " << P(p, 0) << P(p, 1) << P(p, 2) << std::endl; };
		double i = P(p, 0) - corner(0);
		double j = P(p, 1) - corner(1);
		double k = P(p, 2) - corner(2);
		//if (p == 0) { std::cout << "translated point" << i << j << k << std::endl; };

		int ni = i/h;
		int nj = j/h;
		int nk = k/h;
		//std::cout << "cell corner indices: " << ni << nj << nk << std::endl;

		double ci = ni * h;
		double cj = nj * h;
		double ck = nk * h;
		//std::cout << "cell corner position: " << ci << cj << ck << std::endl;

		double offsets[8][3] = { {0, 0, 0},{0, 0, h},{0, h, 0},{0, h, h},{h, 0, 0},{h, 0, h},{h, h, 0},{h, h, h} };
		for (int v = 0; v < 8; v++) {
			double w = (1 - abs(ci + offsets[v][0] - i)/h)* (1 - abs(cj + offsets[v][1] - j)/h)* (1 - abs(ck + offsets[v][2] - k)/h);
			//if (p==0){ std::cout << "W*x: " << w * (ci + offsets[v][0]) << " " << w * (cj + offsets[v][1]) << " " << w * (ck + offsets[v][2]) << std::endl; }
			tripletList.push_back(T(p, (ni + offsets[v][0] / h + (nj + offsets[v][1] / h) * nx + (nk + offsets[v][2] / h) * ny * nx), w));
		}
	}

	W.setFromTriplets(tripletList.begin(), tripletList.end());
}
