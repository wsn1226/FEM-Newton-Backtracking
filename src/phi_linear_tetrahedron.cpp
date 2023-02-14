#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X)
{

    Eigen::Vector3d X0, X1, X2, X3, b, phi3;
    Eigen::MatrixXd T(3, 3);
    X0 << V(element(0), 0), V(element(0), 1), V(element(0), 2);
    X1 << V(element(1), 0), V(element(1), 1), V(element(1), 2);
    X2 << V(element(2), 0), V(element(2), 1), V(element(2), 2);
    X3 << V(element(3), 0), V(element(3), 1), V(element(3), 2);

    b << X - X0;
    T << X1 - X0, X2 - X0, X3 - X0;

    phi3 = T.inverse() * b;
    phi << 1.0 - phi3(0) - phi3(1) - phi3(2), phi3;
}