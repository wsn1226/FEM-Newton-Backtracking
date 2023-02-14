#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>
#include <iostream>
// Input:
//   q - generalized coordinates for the FEM system
//   V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//   element - the 1x4 vertex indices for this tetrahedron
//   v0 - the undeformed tetrahedron volume
//   C,D - material parameters for the Neo-Hookean model
// Output:
//   dV - the 12x12 Hessian of the potential energy for a single tetrahedron
void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q,
                                Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                double C, double D)
{
    Eigen::Matrix3d F;
    Eigen::Matrix34d A;
    Eigen::Matrix43d dphi;
    Eigen::MatrixXd B(9, 12);
    Eigen::Matrix99d d2psi;
    H.setZero();
    B.setZero();
    A.setZero();
    F.setZero();
    dphi.setZero();
    d2psi.setZero();

    Eigen::Vector3d x0, x1, x2, x3, X;
    x0 << q(3 * element(0)), q(3 * element(0) + 1), q(3 * element(0) + 2);
    x1 << q(3 * element(1)), q(3 * element(1) + 1), q(3 * element(1) + 2);
    x2 << q(3 * element(2)), q(3 * element(2) + 1), q(3 * element(2) + 2);
    x3 << q(3 * element(3)), q(3 * element(3) + 1), q(3 * element(3) + 2);

    A << x0, x1, x2, x3;
    dphi_linear_tetrahedron_dX(dphi, V, element, X);
    B.setZero();
    B.block(0, 0, 3, 1) = dphi.block(0, 0, 1, 3).transpose();
    B.block(3, 1, 3, 1) = dphi.block(0, 0, 1, 3).transpose();
    B.block(6, 2, 3, 1) = dphi.block(0, 0, 1, 3).transpose();

    B.block(0, 3, 3, 1) = dphi.block(1, 0, 1, 3).transpose();
    B.block(3, 4, 3, 1) = dphi.block(1, 0, 1, 3).transpose();
    B.block(6, 5, 3, 1) = dphi.block(1, 0, 1, 3).transpose();

    B.block(0, 6, 3, 1) = dphi.block(2, 0, 1, 3).transpose();
    B.block(3, 7, 3, 1) = dphi.block(2, 0, 1, 3).transpose();
    B.block(6, 8, 3, 1) = dphi.block(2, 0, 1, 3).transpose();

    B.block(0, 9, 3, 1) = dphi.block(3, 0, 1, 3).transpose();
    B.block(3, 10, 3, 1) = dphi.block(3, 0, 1, 3).transpose();
    B.block(6, 11, 3, 1) = dphi.block(3, 0, 1, 3).transpose();

    auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X)
    {
        F = A * dphi;
        d2psi_neo_hookean_dF2(d2psi, F, C, D);
        dV = B.transpose() * d2psi * B;
    };
    // integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);

    // DO NOT REMOVE THIS CODE This code ensures that the hessian matrix is symmetric postive definite by projecting all
    // negative eigenvalues to small, postive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);

    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();

    for (int i = 0; i < 12; ++i)
    {
        if (es.eigenvalues()(i) < 1e-6)
        {
            DiagEval(i, i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();
}
