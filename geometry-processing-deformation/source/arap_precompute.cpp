#include "../include/arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/arap_linear_block.h>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>

void arap_precompute(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& b,
    igl::min_quad_with_fixed_data<double>& data,
    Eigen::SparseMatrix<double>& K)
{
    // REPLACE WITH YOUR CODE

    Eigen::SparseMatrix <double> L;
    igl::cotmatrix(V, F, L);

    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute((-L).eval(), b, Aeq, 0, data);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    Eigen::MatrixXd C;
    igl::cotmatrix_entries(V, F, C);
    for (int face = 0; face < F.rows(); face++) {
        for (int edge = 0; edge < 3; edge++) {
            int i = F(face, (edge + 1) % 3);
            int j = F(face, (edge + 2) % 3);
            Eigen::RowVector3d e = C(face, edge) * (V.row(i) - V.row(j));

            for (int ki = 0; ki < 3; ki++) {
                int k = F(face, ki);
                tripletList.push_back(T(i, 3 * k + 0, e(0)));
                tripletList.push_back(T(i, 3 * k + 1, e(1)));
                tripletList.push_back(T(i, 3 * k + 2, e(2)));
                tripletList.push_back(T(j, 3 * k + 0, -e(0)));
                tripletList.push_back(T(j, 3 * k + 1, -e(1)));
                tripletList.push_back(T(j, 3 * k + 2, -e(2)));
            }
        }
    }
    Eigen::SparseMatrix<double> temp(V.rows(), 3 * V.rows());
    temp.setFromTriplets(tripletList.begin(), tripletList.end(), [](const double a, const double b) { return a + b; });

    K = temp/3;

}
