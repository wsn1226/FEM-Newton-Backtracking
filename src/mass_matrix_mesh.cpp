#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

typedef Eigen::Triplet<double> Tri;
void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0)
{
    std::vector<Tri> M_entries;
    Eigen::RowVectorXi ele(4);
    ele.setZero();
    M.resize(qdot.size(), qdot.size());
    M.setZero();

    for (int k = 0; k < T.rows(); k++)
    {
        int ind0 = T(k, 0);
        int ind1 = T(k, 1);
        int ind2 = T(k, 2);
        int ind3 = T(k, 3);
        ele << ind0, ind1, ind2, ind3;
        Eigen::Matrix1212d M0;
        mass_matrix_linear_tetrahedron(M0, qdot, ele, density, v0(k));
        // std::cout << M0(0, 0) << std::endl;
        //-----------------------------FIRST LINE---------------------------------
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                M_entries.push_back(Tri(3 * ind0 + i, 3 * ind0 + j, M0(i, j)));
                M_entries.push_back(Tri(3 * ind0 + i, 3 * ind1 + j, M0(i, j + 3)));
                M_entries.push_back(Tri(3 * ind0 + i, 3 * ind2 + j, M0(i, j + 6)));
                M_entries.push_back(Tri(3 * ind0 + i, 3 * ind3 + j, M0(i, j + 9)));

                //-----------------------------SECOND LINE---------------------------------

                M_entries.push_back(Tri(3 * ind1 + i, 3 * ind0 + j, M0(i + 3, j)));
                M_entries.push_back(Tri(3 * ind1 + i, 3 * ind1 + j, M0(i + 3, j + 3)));
                M_entries.push_back(Tri(3 * ind1 + i, 3 * ind2 + j, M0(i + 3, j + 6)));
                M_entries.push_back(Tri(3 * ind1 + i, 3 * ind3 + j, M0(i + 3, j + 9)));

                //-----------------------------THIRD LINE---------------------------------

                M_entries.push_back(Tri(3 * ind2 + i, 3 * ind0 + j, M0(i + 6, j)));
                M_entries.push_back(Tri(3 * ind2 + i, 3 * ind1 + j, M0(i + 6, j + 3)));
                M_entries.push_back(Tri(3 * ind2 + i, 3 * ind2 + j, M0(i + 6, j + 6)));
                M_entries.push_back(Tri(3 * ind2 + i, 3 * ind3 + j, M0(i + 6, j + 9)));

                //-----------------------------FOURTH LINE---------------------------------

                M_entries.push_back(Tri(3 * ind3 + i, 3 * ind0 + j, M0(i + 9, j)));
                M_entries.push_back(Tri(3 * ind3 + i, 3 * ind1 + j, M0(i + 9, j + 3)));
                M_entries.push_back(Tri(3 * ind3 + i, 3 * ind2 + j, M0(i + 9, j + 6)));
                M_entries.push_back(Tri(3 * ind3 + i, 3 * ind3 + j, M0(i + 9, j + 9)));
            }
        }
    }
    // std::cout << "heyhey" << std::endl;
    M.setFromTriplets(M_entries.begin(), M_entries.end());
}