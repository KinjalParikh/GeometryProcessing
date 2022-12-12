#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
        const Eigen::MatrixXd & l_sqr,
        Eigen::MatrixXd & A)
{
    // Add with your code
    A.resizeLike(l_sqr);
    for (int find = 0; find < l_sqr.rows(); find++) {
        for (int ind = 0; ind < 3; ind++) {
            int ai = (ind + 1) % 3;
            int bi = (ind + 2) % 3;
            int ci = ind;
            double cos = (l_sqr(find, ai) + l_sqr(find, bi) - l_sqr(find, ci)) / (2 * sqrt(l_sqr(find, ai) * l_sqr(find, bi)));
            A(find, ind) = acos(cos);
        }
    }
}