#include "vector_area_matrix.h"
#include "igl/boundary_facets.h"

void vector_area_matrix(
        const Eigen::MatrixXi & F,
        Eigen::SparseMatrix<double>& A)
{
    // Replace with your code
    int V_size = F.maxCoeff()+1;
    A.resize(V_size*2,V_size*2);

    Eigen::MatrixXd bE;
    igl::boundary_facets(F, bE);
    //igl::boundary_facets(F, bE, Eigen::VectorXi(), Eigen::VectorXi());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    Eigen::SparseMatrix<double> A_temp(V_size * 2, V_size * 2);
    for (int eind = 0; eind < bE.rows(); eind++) {
        tripletList.push_back(T(bE(eind, 0), bE(eind, 1) + V_size, 1));
        tripletList.push_back(T(bE(eind, 0) + V_size, bE(eind, 1), -1));
    }

    A_temp.setFromTriplets(tripletList.begin(), tripletList.end());
    A = (A_temp + Eigen::SparseMatrix<double>(A_temp.transpose())) / 2;
}