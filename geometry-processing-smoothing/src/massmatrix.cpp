#include "massmatrix.h"

void massmatrix(
        const Eigen::MatrixXd& l,
        const Eigen::MatrixXi& F,
        Eigen::DiagonalMatrix<double, Eigen::Dynamic>& M)
{
    int num_vertex = F.maxCoeff() + 1;
    Eigen::VectorXd diags(num_vertex);
    diags.setZero();

    for (int f = 0; f < F.rows(); f++) {
        double s = (l(f, 0) + l(f, 1) + l(f, 2)) / 2;
        double area = sqrt(s * (s - l(f, 0)) * (s - l(f, 1)) * (s - l(f, 2)));

        for (int v = 0; v < 3; v++) {
            diags(F(f, v)) += area;
        }
    }

    diags = diags / 3;

    M = diags.asDiagonal();

}