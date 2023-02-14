#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>

// Input:
//   F - the dense 3x3 deformation gradient
//   C,D - material parameters for the Neo-Hookean model
// Output:
//   psi - the neohookean energy
void psi_neo_hookean(double &psi,
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D)
{
    // double J = F.determinant();
    // psi = C * (std::cbrt(std::pow(J, -2)) * ((F.transpose() * F).trace()) - 3.0) + D * std::pow(J - 1.0, 2);
    //   std::cout << "mu:" << mu << std::endl;
    //   std::cout << "lambda:" << lambda << std::endl;
    //    C = mu/2, D = lambda/2
    //   double J = F.determinant();
    //   psi = C * (pow(J, -2.0/3.0) * ((F.transpose()*F).trace()) - 3.0) + D * (J - 1.0)*(J - 1.0);

    // double Ic = F.squaredNorm();
    // double Jminus1 = F.determinant() - 1.0 - mu/lambda;
    // psi = 0.5 * (mu * (Ic - 3.0) + lambda * Jminus1 * Jminus1);
    double mu = 2 * C;
    double lambda = 2 * D;
    double F1_1, F1_2, F1_3, F2_1, F2_2, F2_3, F3_1, F3_2, F3_3;
    F1_1 = F(0, 0);
    F1_2 = F(0, 1);
    F1_3 = F(0, 2);
    F2_1 = F(1, 0);
    F2_2 = F(1, 1);
    F2_3 = F(1, 2);
    F3_1 = F(2, 0);
    F3_2 = F(2, 1);
    F3_3 = F(2, 2);
    psi = (lambda * pow(mu / lambda - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0, 2.0)) / 2.0 + (mu * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3 - 3.0)) / 2.0;
}