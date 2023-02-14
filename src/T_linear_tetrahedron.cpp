#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

// Input:
//  qdot - generalied velocity of FEM system
//  element - vertex indices of the element
//  density - material density
//  volume - volume of tetrahedron
// Output:
//   T - kinetic energy of tetrahedron
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume)
{
    Eigen::Vector12d v;
    Eigen::Matrix1212d M0;
    M0.setZero();

    v << qdot(3 * element(0)), qdot(3 * element(0) + 1), qdot(3 * element(0) + 2),
        qdot(3 * element(1)), qdot(3 * element(1) + 1), qdot(3 * element(1) + 2),
        qdot(3 * element(2)), qdot(3 * element(2) + 1), qdot(3 * element(2) + 2),
        qdot(3 * element(3)), qdot(3 * element(3) + 1), qdot(3 * element(3) + 2);

    mass_matrix_linear_tetrahedron(M0, v, element, density, volume);

    T = 0.5 * v.transpose() * M0 * v;
}