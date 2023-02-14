#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>

// Input:
//   q - generalized coordinates of the FEM system
//   element - vertex indices for the tetrahedron
//  volume - volume of the tetrahedron
//  integrand(out, q, X) - function to be integrated, returns value in out.
// Output:
//   integrated - the value of the integrated function
template <typename Ret, typename Integrand_Func>
// Ret &&integrated
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q,
                                    Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                    Integrand_Func integrand)
{
    Eigen::Vector3d x0, x1, x2, x3, x;
    x0 << q(3 * element(0)), q(3 * element(0) + 1), q(3 * element(0) + 2);
    x1 << q(3 * element(1)), q(3 * element(1) + 1), q(3 * element(1) + 2);
    x2 << q(3 * element(2)), q(3 * element(2) + 1), q(3 * element(2) + 2);
    x3 << q(3 * element(3)), q(3 * element(3) + 1), q(3 * element(3) + 2);

    x = 0.25 * (x0 + x1 + x2 + x3);

    integrand(integrated, q, element, x);
    integrated = integrated * volume;
}
