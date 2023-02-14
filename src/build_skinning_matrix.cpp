#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

// Input:
//   V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//   T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//   V_skin - lx3 matrix of vertices of the display mesh
// Output:
//   N - the lxn sparse skinning matrix

typedef Eigen::Triplet<double> Tri;

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T,
                           Eigen::Ref<const Eigen::MatrixXd> V_skin)
{
    N.resize(V_skin.rows(), V.rows());
    std::vector<Tri> N_entries;
    for (int i = 0; i < V_skin.rows(); i++)
    {
        Eigen::Vector3d X;
        X << V_skin(i, 0), V_skin(i, 1), V_skin(i, 2);
        // run here
        Eigen::RowVectorXi ele(4);
        Eigen::Vector4d phis; // phi0, phi1, phi2, phi3
        // find the tet that contains this vertex, then compute the phis, and put it inside corresponding indices
        for (int j = 0; j < T.rows(); j++)
        {
            int ind0 = T(j, 0);
            int ind1 = T(j, 1);
            int ind2 = T(j, 2);
            int ind3 = T(j, 3);
            ele << ind0, ind1, ind2, ind3;
            // but not run here
            phi_linear_tetrahedron(phis, V, ele, X);
            if (phis(0) > 0 && phis(1) > 0 && phis(2) > 0 && phis(3) > 0)
            {
                N_entries.push_back(Tri(i, ind0, phis(0)));
                N_entries.push_back(Tri(i, ind1, phis(1)));
                N_entries.push_back(Tri(i, ind2, phis(2)));
                N_entries.push_back(Tri(i, ind3, phis(3)));
            }
        }
    }
    N.setFromTriplets(N_entries.begin(), N_entries.end());
}