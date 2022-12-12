#include "closest_rotation.h"
#include "Eigen/dense"

void closest_rotation(
        const Eigen::Matrix3d & M,
        Eigen::Matrix3d & R)
{
    // Replace with your code
    Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeFullU | Eigen::ComputeFullV> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    svd.singularValues();
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();

    double det = (U * V.transpose()).determinant();
    Eigen::Matrix3d omega = Eigen::Matrix3d::Identity();
    omega(2, 2) = det;

    R = U *omega * V.transpose();
}