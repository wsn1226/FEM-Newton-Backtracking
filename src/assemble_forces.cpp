#include <assemble_forces.h>
#include <iostream>

// Input:
//   q - generalized coordinates for the FEM system
//   qdot - generalized velocity for the FEM system
//   V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//   T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//   v0 - the mx1 vector of undeformed tetrahedron volumes
//   C,D - material parameters for the Neo-Hookean model
// Output:
//   f - the vector 3xn vector of forces acting on each node of the mass-spring system

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D)
{
    Eigen::RowVectorXi element(4);
    Eigen::Vector12d dVdq;
    dVdq.setZero();
    f.resize(q.rows());
    f.setZero();
    for (int i = 0; i < T.rows(); i++)
    {
        element << T(i, 0), T(i, 1), T(i, 2), T(i, 3);
        dV_linear_tetrahedron_dq(dVdq, q, V, element, v0(i), C, D);
        f(3 * element(0)) -= dVdq(0);
        f(3 * element(0) + 1) -= dVdq(1);
        f(3 * element(0) + 2) -= dVdq(2);

        f(3 * element(1)) -= dVdq(3);
        f(3 * element(1) + 1) -= dVdq(4);
        f(3 * element(1) + 2) -= dVdq(5);

        f(3 * element(2)) -= dVdq(6);
        f(3 * element(2) + 1) -= dVdq(7);
        f(3 * element(2) + 2) -= dVdq(8);

        f(3 * element(3)) -= dVdq(9);
        f(3 * element(3) + 1) -= dVdq(10);
        f(3 * element(3) + 2) -= dVdq(11);
    }
};