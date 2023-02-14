#include <dV_linear_tetrahedron_dq.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q,
                              Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                              double C, double D)
{
    // first get B (9,12), which needs dphi's entries
    Eigen::Matrix3d F;
    Eigen::Matrix34d A;
    Eigen::Matrix43d dphi;
    Eigen::MatrixXd B(9, 12);
    Eigen::Vector9d dpsi;
    B.setZero();
    A.setZero();
    F.setZero();
    dphi.setZero();

    Eigen::Vector3d x0, x1, x2, x3, X;
    x0 << q(3 * element(0)), q(3 * element(0) + 1), q(3 * element(0) + 2);
    x1 << q(3 * element(1)), q(3 * element(1) + 1), q(3 * element(1) + 2);
    x2 << q(3 * element(2)), q(3 * element(2) + 1), q(3 * element(2) + 2);
    x3 << q(3 * element(3)), q(3 * element(3) + 1), q(3 * element(3) + 2);

    A << x0, x1, x2, x3;

    dphi_linear_tetrahedron_dX(dphi, V, element, X);
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

    // dphi 4,3; F 3,3;
    auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X)
    {
        F = A * dphi;
        dpsi_neo_hookean_dF(dpsi, F, C, D);
        dV = B.transpose() * dpsi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);
    // std::cout << dV << std::endl;
}