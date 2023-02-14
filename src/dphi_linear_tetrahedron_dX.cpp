#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X)
{

    Eigen::Vector3d X0, X1, X2, X3, b, phi3, ONE;
    ONE << 1.0, 1.0, 1.0;
    Eigen::MatrixXd T(3, 3);
    X0 << V(element(0), 0), V(element(0), 1), V(element(0), 2);
    X1 << V(element(1), 0), V(element(1), 1), V(element(1), 2);
    X2 << V(element(2), 0), V(element(2), 1), V(element(2), 2);
    X3 << V(element(3), 0), V(element(3), 1), V(element(3), 2);
    T << X1 - X0, X2 - X0, X3 - X0;
    Eigen::Matrix3d T_inv = T.inverse();
    dphi << -1.0 * ONE.transpose() * T_inv, T_inv;
}