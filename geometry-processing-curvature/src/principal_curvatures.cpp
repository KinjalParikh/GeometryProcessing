#include "../include/principal_curvatures.h"
#include "igl/adjacency_list.h"
#include <set>
#include <igl/eigs.h>
#include <igl/sort.h>
#include <igl/pinv.h>

void principal_curvatures(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        Eigen::MatrixXd & D1,
        Eigen::MatrixXd & D2,
        Eigen::VectorXd & K1,
        Eigen::VectorXd & K2)
{
    // Replace with your code
    K1 = Eigen::VectorXd::Zero(V.rows());
    K2 = Eigen::VectorXd::Zero(V.rows());
    D1 = Eigen::MatrixXd::Zero(V.rows(),3);
    D2 = Eigen::MatrixXd::Zero(V.rows(),3);

    std::vector<std::vector<int>> adj;
    igl::adjacency_list(F, adj);

    for (int vi = 0; vi < V.rows(); vi++) {

        //TWO RING NEIGHBOR
        std::set<int> neighbors;
        std::vector<int> adjvi = adj[vi];
        for (int ai = 0; ai < adjvi.size(); ai++) {
            std::copy(adj[adjvi[ai]].begin(), adj[adjvi[ai]].end(), std::inserter(neighbors, neighbors.end()));
        }

        //RELATIVE POSITION MATRIX P
        Eigen::MatrixXd P(neighbors.size(), 3);
        int index = 0;
        for (int ni : neighbors) {
            P.row(index) = V.row(ni) - V.row(vi);
            index++;
        }

        //COVARIANCE
        Eigen::MatrixXd cov;
        cov = P.transpose() * P;

        //EIGEN DECOMPOSITION TO FIND TANGENT PLANE
        Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
        Eigen::MatrixXd axes = es.eigenvectors().real().eval();
        Eigen::VectorXd vals = es.eigenvalues().real().eval();
        Eigen::VectorXd vals_sorted;
        Eigen::VectorXd map;
        igl::sort(vals, 1, 0, vals_sorted, map);
        Eigen::VectorXd pca1 = axes.col(map[0]);
        Eigen::VectorXd pca2 = axes.col(map[1]);
        Eigen::VectorXd pca3 = axes.col(map[2]);


        //P in u, v coords and w height
        Eigen::MatrixXd UV(neighbors.size(), 2);
        Eigen::VectorXd W(neighbors.size());
        for (int ni = 0; ni<neighbors.size(); ni++) {
            UV(ni, 0) = pca1.dot(P.row(ni));
            UV(ni, 1) = pca2.dot(P.row(ni));
            W(ni) = pca3.dot(P.row(ni));
        }

        //QUADRATIC CURVE FITTING
        Eigen::MatrixXd pinvcoeffs(neighbors.size(), 5);
        Eigen::MatrixXd pinvans;
        pinvcoeffs.col(0) = UV.col(0);
        pinvcoeffs.col(2) = UV.col(0).cwiseProduct(UV.col(0));
        pinvcoeffs.col(1) = UV.col(1);
        pinvcoeffs.col(4) = UV.col(1).cwiseProduct(UV.col(1));
        pinvcoeffs.col(3) = UV.col(0).cwiseProduct(UV.col(1));
        igl::pinv(pinvcoeffs, pinvans);
        Eigen::VectorXd quadcoeffs;
        quadcoeffs = pinvans * W;

        //BUILD SHAPE OPERATOR
        Eigen::Matrix2d S, s1, s2;
        double e, f, g, E, F, G;
        E = 1 + quadcoeffs(0) * quadcoeffs(0);
        F = quadcoeffs(0) * quadcoeffs(1);
        G = 1 + quadcoeffs(1) * quadcoeffs(1);
        e = 2 * quadcoeffs(2) / sqrt(1 + quadcoeffs(0) * quadcoeffs(0) + quadcoeffs(1) * quadcoeffs(1));
        f = quadcoeffs(3) / sqrt(1 + quadcoeffs(0) * quadcoeffs(0) + quadcoeffs(1) * quadcoeffs(1));
        g = 2 * quadcoeffs(4) / sqrt(1 + quadcoeffs(0) * quadcoeffs(0) + quadcoeffs(1) * quadcoeffs(1));
        s1 << e, f, f, g;
        s2 << E, F, F, G;
        S = -s1 * s2.inverse();

        //EIGEN DECOMPOSITION TO FIND PRINCIPLE CURVATURE
        Eigen::EigenSolver<Eigen::MatrixXd> es2(S);
        Eigen::MatrixXd axes2 = es2.eigenvectors().real().eval();
        Eigen::VectorXd vals2 = es2.eigenvalues().real().eval();
        Eigen::VectorXd vals_sorted2;
        Eigen::VectorXd map2;
        igl::sort(vals2, 1, 0, vals_sorted2, map2);
        Eigen::VectorXd pca12 = axes2.col(map2[0]);
        Eigen::VectorXd pca22 = axes2.col(map2[1]);

        //FINAL VALUES
        K1(vi) = vals_sorted2(0);
        K2(vi) = vals_sorted2(1);

        D1.row(vi) = pca12(0) * pca1 + pca12(1) * pca2;
        D2.row(vi) = pca22(0) * pca1 + pca22(1) * pca2;
    }

}