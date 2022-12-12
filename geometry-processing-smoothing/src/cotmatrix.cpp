#include "cotmatrix.h"
#include <iostream>

void cotmatrix(
        const Eigen::MatrixXd& l,
        const Eigen::MatrixXi& F,
        Eigen::SparseMatrix<double>& L)
{

    // Add your code here

    std::map<std::pair<int, int>, double> cots;

    for (int face = 0; face < F.rows(); face++) {

        double s = (l(face, 0) + l(face, 1) + l(face, 2)) / 2;
        double area = sqrt(s * (s - l(face, 0)) * (s - l(face, 1)) * (s - l(face, 2)));

        for (int vs = 0; vs < 3; vs++) {
            int i = F(face, vs);
            int j = F(face, (vs + 1) % 3);
            double la = l(face, (vs + 2) % 3);
            double lb = l(face, (vs + 1) % 3);
            double lc = l(face, vs);
            double cotij = (lb * lb - lc * lc - la * la) / (4 * area);

            if (i > j) {
                int temp = i;
                i = j;
                j = temp;
            }

            auto search = cots.find(std::pair<int, int>(i, j));
            if (search != cots.end()) {
                search->second = search->second + cotij;
            }
            else {
                cots.insert({ std::pair<int, int>(i, j), cotij });
            }

        }
    }

    int num_verts = F.maxCoeff() + 1;
    Eigen::VectorXd diags(num_verts);
    diags.setZero();

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    //int count = 0;
    for (auto const it : cots) {
        double val = it.second / 2;
        if (std::isinf(val)) {
            //std::cout << "\n isinf at index " << it.first.first << " " << it.first.second << std::endl;
            val = std::numeric_limits<double>::max();
            //count += 4;
        }
        else if (std::isnan(val)) {
            //std::cout << "\n isnan at index " << it.first.first << " " << it.first.second << std::endl;
            //count += 4;
            val = 0;
        }

        tripletList.push_back(T(it.first.first, it.first.second, val));
        tripletList.push_back(T(it.first.second, it.first.first, val));
        diags(it.first.first) += val;
        diags(it.first.second) += val;
    }
    //std::cout << "Total degenerate vals " << count << std::endl;

    for (int vert = 0; vert < num_verts; vert++) {
        tripletList.push_back(T(vert, vert, -diags(vert)));
    }

    L.resize(num_verts, num_verts);
    L.setFromTriplets(tripletList.begin(), tripletList.end());

    L = -L;
}