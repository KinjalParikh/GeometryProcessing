#include "tutte.h"
#include "igl/boundary_loop.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/unique_edge_map.h"


void tutte(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        Eigen::MatrixXd & U)
{
    // Replace with your code
    U = V.leftCols(2);
    Eigen::VectorXi bI;
    igl::boundary_loop(F, bI);
    Eigen::MatrixXd bC(bI.size(), 2);
    igl::map_vertices_to_circle(V, bI, bC);

    Eigen::SparseMatrix<double> L(V.rows(), V.rows());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList_un;
    //std::vector<T> tripletList_sum;

    Eigen::MatrixXd uE;
    Eigen::MatrixXd E, emap;
    igl::unique_edge_map(F, E, uE, emap);


    for (int ind = 0; ind < uE.rows(); ind++) {
        int i = uE(ind, 0);
        int j = uE(ind, 1);
        double wij = 1 / (V.row(i) - V.row(j)).norm();

        tripletList_un.push_back(T(i, j, wij));
        tripletList_un.push_back(T(j, i, wij));
        tripletList_un.push_back(T(i, i, -wij));
        tripletList_un.push_back(T(j, j, -wij));
    }

    L.setFromTriplets(tripletList_un.begin(), tripletList_un.end(), [](const double a, const double b) {return a+b; });

    //Eigen::MatrixXd utemp;
    igl::min_quad_with_fixed(L, Eigen::VectorXd::Zero(V.rows()), bI, bC, Eigen::SparseMatrix<double>(), Eigen::MatrixXd(), 0, U);
    U.col(0) = -U.col(0);
}