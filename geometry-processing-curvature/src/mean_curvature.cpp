#include "../include/mean_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/invert_diag.h"
#include "igl/per_vertex_normals.h"
#include "igl/dot_row.h"

void mean_curvature(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        Eigen::VectorXd & H)
{
    // Replace with your code
    H = Eigen::VectorXd::Zero(V.rows());

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    Eigen::MatrixXd Minv;
    igl::invert_diag(M, Minv);

    Eigen::MatrixXd HN;
    HN = - Minv * L * V;

    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);

    H = igl::dot_row(HN, N);
}