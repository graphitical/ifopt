#include <iostream>
#include <iomanip>

#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/test2_vars_constr_cost.h>

using namespace ifopt;

int main()
{
    // 1. Define the problem
    Problem nlp;
    nlp.AddVariableSet  (std::make_shared<ExIpoptVariables>());
    nlp.AddConstraintSet(std::make_shared<ExIpoptConstraint1>());
    nlp.AddConstraintSet(std::make_shared<ExIpoptConstraint2>());
    nlp.AddCostSet      (std::make_shared<ExIpoptCost>());

    // 2. Choose solver and options
    IpoptSolver ipopt;
    ipopt.SetOption("linear_solver","mumps");
    ipopt.SetOption("jacobian_approximation", "exact");
    ipopt.SetOption("hessian_approximation", "exact");
    nlp.PrintCurrent();

    std::cout << std::setprecision(5);
    // 3. Solve
    ipopt.Solve(nlp);
    Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
    std::cout << x.transpose() << std::endl;

    // 4. Test if solution is correct
    double eps = 1e-3; // double precision
    Eigen::Vector4d sol(1.0, 4.74299963, 3.82114998, 1.37940829);
    assert((sol - x).norm() < eps);
}