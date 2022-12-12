#include "../include/angle_defect.h"
#include "igl/squared_edge_lengths.h"
#include "internal_angles.h"

void angle_defect(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        Eigen::VectorXd & D)
{
    D = Eigen::VectorXd::Constant(V.rows(), 2*EIGEN_PI);

    Eigen::MatrixXd L;
    igl::squared_edge_lengths(V, F, L);

    Eigen::MatrixXd A;
    internal_angles(L, A);

    for (int fi = 0; fi < F.rows(); fi++) {
        for (int vi = 0; vi < 3; vi++) {
            D(F(fi, vi)) -= A(fi, vi);
        }
    }

}