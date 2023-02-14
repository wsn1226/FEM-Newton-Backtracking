
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume)
{
    Eigen::Matrix3d I;
    I.setIdentity();
    I /= 120.0;
    M.setZero();

    M << I * 2.0, I, I, I,
        I, I * 2.0, I, I,
        I, I, I * 2.0, I,
        I, I, I, I * 2.0;

    M *= (6.0 * density * volume);
}