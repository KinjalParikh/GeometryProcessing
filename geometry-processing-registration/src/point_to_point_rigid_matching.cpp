#include "point_to_point_rigid_matching.h"
#include "closest_rotation.h"

void point_to_point_rigid_matching(
        const Eigen::MatrixXd & X,
        const Eigen::MatrixXd & P,
        Eigen::Matrix3d & R,
        Eigen::RowVector3d & t)
{
    // Replace with your code
    //R = Eigen::Matrix3d::Identity();
    //t = Eigen::RowVector3d::Zero();


    Eigen::RowVector3d cent_x = X.colwise().sum() / X.rows();
    Eigen::RowVector3d cent_p = P.colwise().sum() / P.rows();

    Eigen::MatrixXd neg_rel_x = X;
    neg_rel_x.rowwise() -= cent_x;
    Eigen::MatrixXd neg_rel_p = P;
    neg_rel_p.rowwise() -= cent_p;

    Eigen::Matrix3d M = neg_rel_p.transpose().eval() * neg_rel_x;
    closest_rotation(M.transpose().eval(), R);

    t = cent_p - (R * cent_x.transpose()).transpose();
}