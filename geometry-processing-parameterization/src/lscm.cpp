#include "lscm.h"
#include "igl/cotmatrix.h"
#include "vector_area_matrix.h"
#include "igl/repdiag.h"
#include "igl/massmatrix.h"
#include "igl/eigs.h"
#include <Eigen/dense>

void lscm(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        Eigen::MatrixXd & U)
{
    // Replace with your code
    U = V.leftCols(2);

    Eigen::SparseMatrix<double> A;
    vector_area_matrix(F, A);

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    Eigen::SparseMatrix<double> Q;
    Q = -igl::repdiag(L, 2) - A;

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);


    Eigen::MatrixXd sU; //#A by k list of sorted eigen vectors
    Eigen::VectorXd sS; //k list of sorted eigen values
    igl::eigs(Q, igl::repdiag(M, 2), 3, igl::EIGS_TYPE_SM, sU, sS);

    U.col(0) = sU.col(2).head(V.rows());
    U.col(1) = sU.col(2).tail(V.rows());

    /*Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeFullV> svd(U, Eigen::ComputeFullV);
    svd.singularValues();
    Eigen::Matrix2d rot = svd.matrixV();

    U = (rot * U.transpose().eval()).transpose().eval();*/


    Eigen::BDCSVD<Eigen::MatrixXd> svd(U.transpose()*U, Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto Urot = svd.matrixU();
    auto Vrot = svd.matrixV();
    U = U * Vrot;

}