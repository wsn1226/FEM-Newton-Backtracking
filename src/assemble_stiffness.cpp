#include <assemble_stiffness.h>
#include <iostream>

typedef Eigen::Triplet<double> Tri;
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                        double C, double D)
{
    std::vector<Tri> K_entries;
    Eigen::RowVectorXi ele(4);
    ele.setZero();
    K.resize(qdot.size(), qdot.size());
    K.setZero();

    for (int k = 0; k < T.rows(); k++)
    {
        int ind0 = T(k, 0);
        int ind1 = T(k, 1);
        int ind2 = T(k, 2);
        int ind3 = T(k, 3);
        ele << ind0, ind1, ind2, ind3;
        Eigen::Matrix1212d K0;
        d2V_linear_tetrahedron_dq2(K0, q, V, ele, v0(k), C, D);
        // std::cout << K0.coeff(0, 0) << std::endl;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                //-----------------------------FIRST LINE---------------------------------
                K_entries.push_back(Tri(3 * ind0 + i, 3 * ind0 + j, K0(i, j)));
                K_entries.push_back(Tri(3 * ind0 + i, 3 * ind1 + j, K0(i, j + 3)));
                K_entries.push_back(Tri(3 * ind0 + i, 3 * ind2 + j, K0(i, j + 6)));
                K_entries.push_back(Tri(3 * ind0 + i, 3 * ind3 + j, K0(i, j + 9)));

                //-----------------------------SECOND LINE---------------------------------

                K_entries.push_back(Tri(3 * ind1 + i, 3 * ind0 + j, K0(i + 3, j)));
                K_entries.push_back(Tri(3 * ind1 + i, 3 * ind1 + j, K0(i + 3, j + 3)));
                K_entries.push_back(Tri(3 * ind1 + i, 3 * ind2 + j, K0(i + 3, j + 6)));
                K_entries.push_back(Tri(3 * ind1 + i, 3 * ind3 + j, K0(i + 3, j + 9)));

                //-----------------------------THIRD LINE---------------------------------

                K_entries.push_back(Tri(3 * ind2 + i, 3 * ind0 + j, K0(i + 6, j)));
                K_entries.push_back(Tri(3 * ind2 + i, 3 * ind1 + j, K0(i + 6, j + 3)));
                K_entries.push_back(Tri(3 * ind2 + i, 3 * ind2 + j, K0(i + 6, j + 6)));
                K_entries.push_back(Tri(3 * ind2 + i, 3 * ind3 + j, K0(i + 6, j + 9)));

                //-----------------------------FOURTH LINE---------------------------------

                K_entries.push_back(Tri(3 * ind3 + i, 3 * ind0 + j, K0(i + 9, j)));
                K_entries.push_back(Tri(3 * ind3 + i, 3 * ind1 + j, K0(i + 9, j + 3)));
                K_entries.push_back(Tri(3 * ind3 + i, 3 * ind2 + j, K0(i + 9, j + 6)));
                K_entries.push_back(Tri(3 * ind3 + i, 3 * ind3 + j, K0(i + 9, j + 9)));
            }
        }
    }
    // std::cout << "heyhey" << std::endl;
    K.setFromTriplets(K_entries.begin(), K_entries.end());
    K *= -1.0;
};
