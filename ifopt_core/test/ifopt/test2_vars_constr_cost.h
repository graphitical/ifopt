
/**
 *  @file test2_vars_constr_cost.h
 *
 *  @brief Example to generate a solver-independent formulation for the problem, taken
 *  from the IPOPT interfaces example: https://coin-or.github.io/Ipopt/INTERFACES.html
 *  
 *
 *  The example problem to be solved is given as:
 *
 *      min_x f(x) = x1*x4 * (x1 + x2 + x3) + x3
 *      s.t.
 *          x1*x2*x3*x4 >= 25
 *          x1^2 + x2^2 + x3^2 + x4^2 = 40
 *          1 <= x1, x2, x3, x4 <= 5
 * 
 *          x_start = (1, 5, 5, 1)
 *          x* = (1.0000, 4.74300, 3.82115, 1.37941)
 * 
 * In this simple example we use one set of variables, two constraints and one
 * cost; however, most real world problems have multiple different constraints
 * and also different variable sets representing different quantities. This
 * framework allows to define each set of variables or constraints absolutely
 * independently from another and correctly stitches them together to form the
 * final optimization problem.
 *
 * For a helpful graphical overview, see:
 * http://docs.ros.org/api/ifopt/html/group__ProblemFormulation.html
 */


#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>

namespace ifopt {
using Eigen::Vector4d;

class ExIpoptVariables : public VariableSet {
public:
    ExIpoptVariables() : ExIpoptVariables("var_set1") {};
    ExIpoptVariables(const std::string& name) : VariableSet(4, name)
    {
        x_(0) = 1.;
        x_(1) = 5.;
        x_(2) = 5.;
        x_(3) = 1.;
    }

    void SetVariables(const VectorXd& x) override
    {
        assert(x.size() == 4);
        x_(0) = x(0);
        x_(1) = x(1);
        x_(2) = x(2);
        x_(3) = x(3);
    }

    VectorXd GetValues() const override
    {
        return Vector4d(x_(0), x_(1), x_(2), x_(3));
    }

    VecBound GetBounds() const override
    {
        VecBound bounds(GetRows());
        for (int i = 0; i < GetRows(); ++i)
            bounds.at(i) = Bounds(1., 5.);
        return bounds;
    }

private:
    Vector4d x_;
};


class ExIpoptConstraint1 : public ConstraintSet {
public:
    ExIpoptConstraint1() : ExIpoptConstraint1("constraint1") {};

    ExIpoptConstraint1(const std::string& name) : ConstraintSet(1, name) {};

    VectorXd GetValues() const override
    {
        VectorXd g(GetRows());
        Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();
        g(0) = x.prod();
        return g;
    }

    VecBound GetBounds() const override
    {
        VecBound b(GetRows());
        b.at(0) = Bounds(25., inf);
        return b;
    }

    void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
    {
        if (var_set == "var_set1") {
            Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();
            int n = x.size();

            for (int i = 0; i < x.size(); ++i)
                jac_block.coeffRef(0, i) = x((i+1)%n) * x((i+2)%n) * x((i+3)%n);
        }
    }

    void FillHessianBlock (std::string var_set, Hessian& hes_block) const override
    {
        if (var_set == "var_set1") {
            Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();
            int n = x.size();

            hes_block.coeffRef(1,0) = x(2) * x(3);
            hes_block.coeffRef(2,0) = x(1) * x(3);
            hes_block.coeffRef(2,1) = x(0) * x(3);
            hes_block.coeffRef(3,0) = x(1) * x(2);
            hes_block.coeffRef(3,1) = x(0) * x(2);
            hes_block.coeffRef(3,2) = x(0) * x(1);
        }
        
    }
};


class ExIpoptConstraint2 : public ConstraintSet {
public:
    ExIpoptConstraint2() : ExIpoptConstraint2("constraint2") {};

    ExIpoptConstraint2(const std::string& name) : ConstraintSet(1, name) {}

    VectorXd GetValues() const override
    {
        VectorXd g(GetRows());
        Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();
        g(0) = x.squaredNorm();
        return g;
    }

    VecBound GetBounds() const override
    {
        VecBound b(GetRows());
        b.at(0) = Bounds(40., 40.);
        return b;
    }

    void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override 
    {
        if (var_set == "var_set1") {
            Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();

            for (int i = 0; i < x.size(); ++i)
                jac_block.coeffRef(0, i) = 2. * x(i);
        }
    }

    void FillHessianBlock (std::string var_set, Hessian& hes_block) const override
    {
        if (var_set == "var_set1") {
            Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();

            for (int i = 0; i < x.size(); ++i)
                hes_block.coeffRef(i,i) = 2.;
        }
    }
};


class ExIpoptCost : public CostTerm {
public: 
    ExIpoptCost() : ExIpoptCost("cost_term1") {};
    ExIpoptCost(const std::string& name) : CostTerm(name) {};

    double GetCost() const override 
    {
        Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();
        return x(0)*x(3) * (x(0) + x(1) + x(2)) + x(2);
    }

    void FillJacobianBlock (std::string var_set, Jacobian& jac) const override
    {
        if (var_set == "var_set1") {
            Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();

            jac.coeffRef(0,0) = x(0)*x(3) + x(3)*(x(0) + x(1) + x(2));
            jac.coeffRef(0,1) = x(0)*x(3);
            jac.coeffRef(0,2) = x(0)*x(3) + 1.;
            jac.coeffRef(0,3) = x(0) * (x(0) + x(1) + x(2));
        }
    }
    void FillHessianBlock (std::string var_set, Hessian& hes) const override
    {
        if (var_set == "var_set1") {
            Vector4d x = GetVariables()->GetComponent("var_set1")->GetValues();

            hes.coeffRef(0,0) = 2. * x(3);
            hes.coeffRef(1,0) = x(3);
            hes.coeffRef(2,0) = x(3);
            hes.coeffRef(3,0) = 2.*x(0) + x(1) + x(2);
            hes.coeffRef(3,1) = x(0);
            hes.coeffRef(3,2) = x(0);
        }
    }

};

} // namespace ifopt