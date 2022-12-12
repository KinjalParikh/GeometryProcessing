#include "point_to_plane_rigid_matching.h"
#include <Eigen/Dense>
#include "closest_rotation.h"

void point_to_plane_rigid_matching(
        const Eigen::MatrixXd & X,
        const Eigen::MatrixXd & P,
        const Eigen::MatrixXd & N,
        Eigen::Matrix3d & R,
        Eigen::RowVector3d & t)
{
    // Replace with your code

    Eigen::MatrixXd A(6, 6);
    A.setZero();
    Eigen::VectorXd b(6);
    b.setZero();


    for (int i = 0; i < X.rows(); i++) {
        Eigen::VectorXd xcn;
        Eigen::Vector3d x = X.row(i);
        Eigen::Vector3d n = N.row(i);
        Eigen::Vector3d p = P.row(i);
        xcn = x.cross(n);

        Eigen::VectorXd a1(6);
        a1.head(3) = xcn;
        a1.tail(3) = n;
        Eigen::MatrixXd a2(6, 6);
        a2 = a1 * a1.transpose();
        A += a2;
        b += a1 * n.transpose() * (p - x);
    }

    Eigen::VectorXd u;
    u = A.inverse() * b;

    t << u(3), u(4), u(5);
    Eigen::RowVector3d w;
    w << u(0), u(1), u(2);

    //double theta = aa.norm();
    //Eigen::RowVector3d w = aa/theta;

    Eigen::Matrix3d Rm = Eigen::Matrix3d::Identity();
    Eigen::MatrixXd W(3, 3);
    W << 0, -w(2), w(1),
            w(2), 0, -w(0),
            -w(1), w(0), 0;
    //R += sin(theta) * W + (1 - cos(theta)) * W * W;
    Rm += W;
    closest_rotation(Rm.transpose().eval(), R);
}