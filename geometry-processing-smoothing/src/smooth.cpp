#include "smooth.h"
#include <igl/edge_lengths.h>
#include "massmatrix.h"
#include "cotmatrix.h"
#include <Eigen/SparseCholesky>
#include <iostream>

void smooth(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& G,
        double lambda,
        Eigen::MatrixXd& U)
{
    // Replace with your code
    Eigen::MatrixXd l(F.rows(), 3);
    igl::edge_lengths(V, F, l);

    const int num_v = V.rows();
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
    massmatrix(l, F, M);

    Eigen::SparseMatrix<double> mass2 = Eigen::SparseMatrix<double>(M);
    //std::cout << "equality: " << mass.isApprox(mass2) << std::endl;

    Eigen::SparseMatrix<double> L2;
    cotmatrix(l, F, L2);

    Eigen::MatrixXd b = mass2 * G;
    Eigen::SparseMatrix<double> A = mass2 - lambda * L2;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    U = solver.solve(b);
    //U = G;
}