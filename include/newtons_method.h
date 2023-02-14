#include <Eigen/Dense>
#include <EigenTypes.h>
#include <Eigen/SparseCholesky>

// Input:
//   x0 - initial point for newtons search
//   f(x) - function that evaluates and returns the cost function at x
//   g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//   H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//   max steps - the maximum newton iterations to take
//   tmp_g and tmp_H are scratch space to store gradients and hessians
// Output:
//   x0 - update x0 to new value
template <typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H)
{
   double p, c, alpha, E, tol;
   int ite_num = 0;
   Eigen::VectorXd d;
   Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
   p = 0.5;
   c = 1e-8;
   tol = 1e-4;
   g(tmp_g, x0);
   H(tmp_H, x0);
   while (tmp_g.norm() > tol && ite_num < maxSteps - 1)
   {
      solver.compute(tmp_H);
      d = solver.solve(-1.0 * tmp_g);
      // Backtracking line search to choose alpha, if without we will do x0+=d
      alpha = 1.0;
      E = f(x0 + alpha * d);
      while (alpha > tol && E > f(x0) + alpha * c * d.transpose() * tmp_g)
      {
         alpha *= p;
         E = f(x0 + alpha * d);
      }
      x0 += alpha * d;
      g(tmp_g, x0);
      H(tmp_H, x0);
      ite_num += 1;
   }

   return 0.0;
}
