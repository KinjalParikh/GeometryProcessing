#include "random_points_on_mesh.h"
#include <igl/doublearea.h>
#include <igl/cumsum.h>
#include <random>
#include <iostream>


void random_points_on_mesh(
        const int n,
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        Eigen::MatrixXd & X)
{
    // REPLACE WITH YOUR CODE:
    X.resize(n,3);

    Eigen::VectorXd area;
    igl::doublearea(V, F, area);
    area = 0.5 * area;
    Eigen::VectorXd ac;
    igl::cumsum(area, 1, ac);
    ac = ac / area.sum();

    for (int i = 0; i < n; i++) {

        //find random number between 0 and 1. Github calls this a2
        double a = (double) rand() / RAND_MAX;
        //std::cout << "a2: " << a << std::endl;

        //find the face
        int l = 0;
        int u = F.rows();
        int index = u/2;
        while (true) {
            if (ac(index) == a) {
                if (index + 1 < F.rows()) {
                    index = index + 1;
                }
                break;
            }
            else if (ac(index) > a) {
                if (index==0 || ac(index - 1) < a) break;
                else {
                    u = index;
                }
            }
            else {
                l = index;
            }
            int next_ind = (l + u) / 2;
            if (next_ind == index) index += 1;
            else index = next_ind;
        }
        //std::cout << "index: " << index << std::endl;

        //find the point
        double a1 = (double) rand() / RAND_MAX;
        if (a + a1 > 1) {
            a = 1 - a;
            a1 = 1 - a1;
        }
        //std::cout << "a1: " << a1 << std::endl;
        X.row(i) = V.row(F(index, 0)) + a1 * (V.row(F(index, 1)) - V.row(F(index, 0))) + a * (V.row(F(index, 2)) - V.row(F(index, 0)));
        //std::cout << "point: " << X.row(i) << std::endl;

    }

}