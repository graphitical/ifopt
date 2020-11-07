/******************************************************************************
Copyright (c) 2017, Alexander W Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <ifopt/problem.h>
#include <iostream>
#include <iomanip>


namespace ifopt {

Problem::Problem ()
    :constraints_("constraint-sets", false),
     costs_("cost-terms", true)
{
  variables_ = std::make_shared<Composite>("variable-sets", false);
}

void
Problem::AddVariableSet(VariableSet::Ptr variable_set)
{
  variables_->AddComponent(variable_set);
}

void
Problem::AddConstraintSet(ConstraintSet::Ptr constraint_set)
{
  constraint_set->LinkWithVariables(variables_);
  constraints_.AddComponent(constraint_set);
}

void
Problem::AddCostSet(CostTerm::Ptr cost_set)
{
  cost_set->LinkWithVariables(variables_);
  costs_.AddComponent(cost_set);
}

int
Problem::GetNumberOfOptimizationVariables () const
{
  return variables_->GetRows();
}

Problem::VecBound
Problem::GetBoundsOnOptimizationVariables () const
{
  return variables_->GetBounds();
}

Problem::VectorXd
Problem::GetVariableValues () const
{
  return variables_->GetValues();
}

void
Problem::SetVariables (const double* x)
{
  variables_->SetVariables(ConvertToEigen(x));
}

double
Problem::EvaluateCostFunction (const double* x)
{
  VectorXd g = VectorXd::Zero(1);
  if (HasCostTerms()) {
    SetVariables(x);
    g = costs_.GetValues();
  }
  return g(0);
}

Problem::VectorXd
Problem::EvaluateCostFunctionGradient (const double* x, bool use_finite_difference_approximation)
{
  int n = GetNumberOfOptimizationVariables();
  Jacobian jac = Jacobian(1,n);
  if (HasCostTerms()) {
    if(use_finite_difference_approximation) {
      // ipopt uses 10e-8 for their derivative check.
      // https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_derivative_test_perturbation
      // setting here to more precise, but can be adapted
      double step_size = std::numeric_limits<double>::epsilon();
      
      // calculate forward difference by disturbing each optimization variable
      double g = EvaluateCostFunction(x);
      std::vector<double> x_new(x, x + n);
      for (int i=0; i<n; ++i) {
        x_new[i] += step_size; // disturb
        double g_new = EvaluateCostFunction(x_new.data());
        jac.coeffRef(0,i) = (g_new - g)/step_size;
        x_new[i] -= step_size; // reset for next iteration
      }
    } else {
      SetVariables(x);
      jac = costs_.GetJacobian();
    }
  }

  return jac.row(0).transpose();
}

Problem::Hessian
Problem::EvaluateCostFunctionHessian (const double* x)
{
  int total_n_var = GetNumberOfOptimizationVariables();
  Hessian hes = Hessian(total_n_var, total_n_var);
  if (HasCostTerms()) {
    SetVariables(x);
    hes = costs_.GetHessian();
  }

  // TODO: Does this need to be transposed since we're storing only the bottom left triangle?
  // Hessian assumed to be symmetric so IPOPT only wants the bottom left triangle
  return hes;
}

Problem::VecBound
Problem::GetBoundsOnConstraints () const
{
  return constraints_.GetBounds();
}

int
Problem::GetNumberOfConstraints () const
{
  return GetBoundsOnConstraints().size();
}

Problem::VectorXd
Problem::EvaluateConstraints (const double* x)
{
  SetVariables(x);
  return constraints_.GetValues();
}

bool
Problem::HasCostTerms () const
{
  return costs_.GetRows()>0;
}

void
Problem::EvalNonzerosOfJacobian (const double* x, double* values)
{
  SetVariables(x);
  Jacobian jac = GetJacobianOfConstraints();

  jac.makeCompressed(); // so the valuePtr() is dense and accurate
  std::copy(jac.valuePtr(), jac.valuePtr() + jac.nonZeros(), values);
}

void
Problem::EvalNonzerosOfHessian (const double* x, const double obj_factor,
                                const bool new_lambda, const double* lambda,
                                double* values)
{
  SetVariables(x);
  Hessian cost_hes = GetHessianOfCosts(obj_factor);
  // TODO: Do something with new_lambda so we don't have to recalculate constraints every time
  Hessian cons_hes = GetHessianOfConstraints(lambda);
  std::cout << "EVALUATING HESSIAN NONZEROS\n";
  
  // Debug. TODO: Delete
  // The Hessians should both be n_var by n_var and x should be size n_var
  int n_var = GetNumberOfOptimizationVariables();
  int n_con = GetNumberOfConstraints();
  std::cout << "Cost Hes:\n" << cost_hes.toDense() << std::endl;

  for (int i = 0; i < n_var; ++i)
    std::cout << "var" << i << ": " << x[i] << std::endl;

  std::cout << "Cons Hes:\n" << cons_hes.toDense() << std::endl
            << "new_lambda: " << new_lambda << "\n";
  for (int i = 0; i < n_con; ++i)
    std::cout << "lam" << i << ": " << lambda[i] << std::endl;

  assert(n_var == cost_hes.cols() && n_var == cons_hes.cols());

  Hessian hes(n_var, n_var);
  std::vector< Eigen::Triplet<double> > triplet_list;
  triplet_list.reserve(cost_hes.nonZeros() + cons_hes.nonZeros());

  for (int k=0; k<cost_hes.outerSize(); ++k)
    for (Hessian::InnerIterator it(cost_hes,k); it; ++it)
      triplet_list.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));

  for (int k=0; k<cons_hes.outerSize(); ++k)
    for (Hessian::InnerIterator it(cons_hes,k); it; ++it)
      triplet_list.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));

  hes.setFromTriplets(triplet_list.begin(), triplet_list.end());

  std::cout << "Final Hessian:\n" << hes.toDense() << std::endl;

  hes.makeCompressed(); // so the valuePtr() is dense and accurate
  std::copy(hes.valuePtr(), hes.valuePtr() + hes.nonZeros(), values);
}

Problem::Jacobian
Problem::GetJacobianOfConstraints () const
{
  return constraints_.GetJacobian();
}

Problem::Jacobian
Problem::GetJacobianOfCosts () const
{
  return costs_.GetJacobian();
}

// ??? What happens on first pass? What are default lambda values?
// ??? Instead of returning Hessian, return triplet_list so I can build it all in the EvalNonZerosHessian
Problem::Hessian
Problem::GetHessianOfConstraints (const double* lambda) const
{
  std::cout << "GET CONSTRAINT HESSIAN\n";

  int n_var = constraints_.GetComponents().empty() ? 0 : constraints_.GetComponents().front()->GetHessian().cols();
  Hessian hessian(n_var, n_var);

  if (n_var == 0) return hessian;

  int con = 0;
  double lam = 0.;
  std::vector< Eigen::Triplet<double> > triplet_list;

  for (const auto& c : constraints_.GetComponents()) {
    const Hessian& hes = c->GetHessian();
    triplet_list.reserve(triplet_list.size() + hes.nonZeros());

      if (lambda == nullptr)
        lam = 1.;
      else
        lam = lambda[con];

    for (int k=0; k<hes.outerSize(); ++k)
      for (Hessian::InnerIterator it(hes,k); it; ++it)
        triplet_list.push_back(Eigen::Triplet<double>(it.row(), it.col(), lam * it.value()));

    // ??? Do I need to bump rows for multiple cost function?
    con++;
  }
  hessian.setFromTriplets(triplet_list.begin(), triplet_list.end());
  std::cout << "DONE WITH THE CONSTRAINT HESSIAN\n";
  return hessian;
}

Problem::Hessian
Problem::GetHessianOfCosts (const double obj_factor) const
{
  std::cout << "GET COST HESSIAN\n";
  // return obj_factor * costs_.GetHessian();
  // TODO: Change back to line above
  Hessian hes = obj_factor * costs_.GetHessian();
  std::cout << "DONE WITH COST HESSIAN\n";
  return hes;
}

void
Problem::SaveCurrent()
{
  x_prev.push_back(variables_->GetValues());
}

Composite::Ptr
Problem::GetOptVariables () const
{
  return variables_;
}

void
Problem::SetOptVariables (int iter)
{
  variables_->SetVariables(x_prev.at(iter));
}

void
Problem::SetOptVariablesFinal ()
{
  variables_->SetVariables(x_prev.at(GetIterationCount()-1));
}

void
Problem::PrintCurrent() const
{
  using namespace std;
  cout << "\n"
       << "************************************************************\n"
       << "    IFOPT - Interface to Nonlinear Optimizers (v2.0)\n"
       << "                \u00a9 Alexander W. Winkler\n"
       << "           https://github.com/ethz-adrl/ifopt\n"
       << "************************************************************"
       << "\n"
       << "Legend:\n"
       << "c - number of variables, constraints or cost terms" << std::endl
       << "i - indices of this set in overall problem" << std::endl
       << "v - number of [violated variable- or constraint-bounds] or [cost term value]"
       << "\n\n"
       << std::right
       << std::setw(33) << ""
       << std::setw(5)  << "c  "
       << std::setw(16) << "i    "
       << std::setw(11) << "v "
       << std::left
       << "\n";

  variables_->PrintAll();
  constraints_.PrintAll();
  costs_.PrintAll();
};

Problem::VectorXd
Problem::ConvertToEigen(const double* x) const
{
  return Eigen::Map<const VectorXd>(x,GetNumberOfOptimizationVariables());
}

} /* namespace opt */

