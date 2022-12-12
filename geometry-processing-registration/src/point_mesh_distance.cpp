#include "point_mesh_distance.h"
#include "point_triangle_distance.h"
#include <Eigen/Dense>

void point_mesh_distance(
        const Eigen::MatrixXd & X,
        const Eigen::MatrixXd & VY,
        const Eigen::MatrixXi & FY,
        Eigen::VectorXd & D,
        Eigen::MatrixXd & P,
        Eigen::MatrixXd & N)
{
    // Replace with your code
    D.resize(X.rows());
    P.resizeLike(X);
    N = Eigen::MatrixXd::Zero(X.rows() ,X.cols());

    for (int i = 0; i < X.rows(); i++) {
        double min_dist;
        Eigen::RowVector3d min_point;
        Eigen::RowVector3d min_normal;
        for (int j = 0; j < FY.rows(); j++) {
            double dist;
            Eigen::RowVector3d point;

            Eigen::RowVector3d x = X.row(i);
            Eigen::RowVector3d a = VY.row(FY(j, 0));
            Eigen::RowVector3d b = VY.row(FY(j, 1));
            Eigen::RowVector3d c = VY.row(FY(j, 2));

            point_triangle_distance(x, a, b, c, dist, point);
            if (j == 0 || dist < min_dist) {
                min_dist = dist;
                min_point = point;
                min_normal = (b - a).cross(c - a);
            }
        }
        D(i) = min_dist;
        P.row(i) = min_point;
        N.row(i) = min_normal;

    }

}