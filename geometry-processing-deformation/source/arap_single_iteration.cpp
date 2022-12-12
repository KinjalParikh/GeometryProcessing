#include "../include/arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/fit_rotations.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
	const igl::min_quad_with_fixed_data<double>& data,
	const Eigen::SparseMatrix<double>& K,
	const Eigen::MatrixXd& bc,
	Eigen::MatrixXd& U)
{
	// REPLACE WITH YOUR CODE
	Eigen::MatrixXd C;
	C = (U.transpose() * K).transpose();
	Eigen::MatrixXd R;
	R.resizeLike(C);
	for (int i = 0; i < U.rows(); i++) {
		Eigen::Matrix3d R2;
		Eigen::Matrix3d c2;
		c2 << C.block(i * 3, 0, 3, 3);
		igl::polar_svd3x3(c2, R2);
		R.block(i * 3, 0, 3, 3) = R2;
	}

	Eigen::MatrixXd b = -(R.transpose() * K.transpose()).transpose();
	igl::min_quad_with_fixed_solve(data, b, bc, Eigen::MatrixXd(), U);
}
