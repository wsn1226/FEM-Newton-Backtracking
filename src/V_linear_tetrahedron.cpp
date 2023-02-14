#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

// Input:
//  q - generalized coordinates of FEM system
//  V - vertex matrix for the mesh
//  element - vertex indices of the element
//  volume - volume of tetrahedron
//  C,D - material parameters
// Output:
//   energy - potential energy of tetrahedron

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q,
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D)
{
    Eigen::Matrix3d F;
    Eigen::Matrix34d A;
    Eigen::Matrix43d dphi;
    A.setZero();
    F.setZero();
    dphi.setZero();

    Eigen::Vector3d x0, x1, x2, x3, X0, X1, X2, X3, X;
    x0 << q(3 * element(0)), q(3 * element(0) + 1), q(3 * element(0) + 2);
    x1 << q(3 * element(1)), q(3 * element(1) + 1), q(3 * element(1) + 2);
    x2 << q(3 * element(2)), q(3 * element(2) + 1), q(3 * element(2) + 2);
    x3 << q(3 * element(3)), q(3 * element(3) + 1), q(3 * element(3) + 2);

    A << x0, x1, x2, x3;

    dphi_linear_tetrahedron_dX(dphi, V, element, X);

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X)
    {
        // We need F, F = (x0, x1, x2, x3)*dphi
        F = A * dphi;
        psi_neo_hookean(e, F, C, D);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);
}