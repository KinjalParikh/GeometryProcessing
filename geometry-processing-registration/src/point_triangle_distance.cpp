#include "point_triangle_distance.h"
#include "igl/barycentric_coordinates.h"

void point_triangle_distance(
        const Eigen::RowVector3d & x,
        const Eigen::RowVector3d & a,
        const Eigen::RowVector3d & b,
        const Eigen::RowVector3d & c,
        double & d,
        Eigen::RowVector3d & p)
{
    // Replace with your code
    //d = 0;
    //p = a;

    Eigen::RowVector3d n = (b - a).cross(c-a);
    n.normalize();

    Eigen::RowVector3d point = x + (n.dot(a) - n.dot(x))*(n);
    Eigen::RowVector3d barycords;
    igl::barycentric_coordinates(x, a, b, c, barycords);

    if (barycords.minCoeff() > 0 && barycords.maxCoeff() < 1) {
        p = point;
        d = (p - x).norm();
    }
    else {
        for (int i = 0; i < 3; i++) {
            if (barycords(i) < 0) barycords(i) = 0;
            else if (barycords(i) > 1) barycords(i) = 1;
        }

        barycords /= barycords.sum();

        p = barycords(0) * a + barycords(1) * b + barycords(2) * c;
        d = (p - x).norm();
    }
}