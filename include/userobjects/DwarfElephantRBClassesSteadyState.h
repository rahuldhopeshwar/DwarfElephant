/**
 * In this class subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
 *
 * DwarfElephantRBEvaluation: requires only the definition of the lower
 * coercivity constant. The value is here specified for a three layer
 * problem.
 *
 * DwarfElephantRBConstruction: In order to construct the RB System with the
 * DwarfElephantRBEvaluation subclass the method build_rb_evaluation needs to be
 * overriden.
 *
 * NOTE: ENSURE THAT THE CLASS IS USING THE CORRECT RBSTRUCTURES.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBCLASSESSTEADYSTATE_H
#define DWARFELEPHANTRBCLASSESSTEADYSTATE_H

///---------------------------------INCLUDES--------------------------------
//#if defined(LIBMESH_HAVE_SLEPC) && defined(LIBMESH_HAVE_GLPK)

// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"

// libMesh includes (RB package)
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"
#include "libmesh/rb_scm_construction.h"
#include "libmesh/rb_scm_evaluation.h"

///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "CacheBoundaries.h"
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "RBStructuresP1T1EqualF1O1SteadyState.h"
#include "RBStructuresP1T2EqualF1O1SteadyState.h"
#include "RBStructuresP1T3EqualF1O1SteadyState.h"

#include "FEProblemBase.h"

// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;

  class EquationSystems;
  class RBConstruction;
  class RBEvaluation;
}

//using libMesh::RBSCMConstruction;
//using libMesh::RBSCMEvaluation;

class DwarfElephantInitializeRBSystemSteadyState;

///-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
class DwarfElephantRBConstructionSteadyState : public RBConstruction
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBConstructionSteadyState (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : Parent(es, name_in, number_in)
  {}

  // Destructor
  virtual ~DwarfElephantRBConstructionSteadyState () { }

  // Type of the system
  typedef DwarfElephantRBConstructionSteadyState _sys_type;

  // Type of the parent
  typedef RBConstruction Parent;

  // Initialize data structure
  virtual void init_data()
  {
    u_var = this->add_variable ("u", libMesh::FIRST);

    Parent::init_data();
  }

  Real compute_residual_dual_norm(const unsigned int N)
{
   LOG_SCOPE("compute_residual_dual_norm()", "RBConstruction");

   // Put the residual in rhs in order to compute the norm of the Riesz representor
   // Note that this only works in serial since otherwise each processor will
   // have a different parameter value during the Greedy training.

   UniquePtr< NumericVector<Number> > RB_sol = NumericVector<Number>::build(this->comm());
   RB_sol->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

   UniquePtr< NumericVector<Number> > temp = NumericVector<Number>::build(this->comm());
   temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

   for(unsigned int i=0; i<N; i++)
   {
     RB_sol->add(get_rb_evaluation().RB_solution(i), get_rb_evaluation().get_basis_function(i));
   }

   this->truth_assembly();
   matrix->vector_mult(*temp, *RB_sol);
   rhs->add(-1., *temp);

   // Then solve to get the Reisz representor
   matrix->zero();
   matrix->add(1., *inner_product_matrix);
//   if(constrained_problem)
//     matrix->add(1., *constraint_matrix);

   solution->zero();
   solve();
   // Make sure we didn't max out the number of iterations
   if( (this->n_linear_iterations() >=
        this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
       (this->final_linear_residual() >
        this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
   {
     libmesh_error_msg("Warning: Linear solver may not have converged! Final linear residual = "
                       << this->final_linear_residual() << ", number of iterations = "
                       << this->n_linear_iterations());
   }

   inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);

   Real slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);

   return std::sqrt( libmesh_real(slow_residual_norm_sq) );
}

//Real truth_solve(int plot_solution)
//{
//  LOG_SCOPE("truth_solve()", "RBConstruction");
//
//  truth_assembly();
//
//  // truth_assembly assembles into matrix and rhs, so use those for the solve
//  if (extra_linear_solver)
//    {
//      // If extra_linear_solver has been initialized, then we use it for the
//      // truth solves.
//      solve_for_matrix_and_rhs(*extra_linear_solver, *matrix, *rhs);
//
//      if (assert_convergence)
//        check_convergence(*extra_linear_solver);
//    }
//  else
//    {
//      solve_for_matrix_and_rhs(*get_linear_solver(), *matrix, *rhs);
//
//      if (assert_convergence)
//        check_convergence(*get_linear_solver());
//    }
//
//
//
//  const RBParameters & mu = get_parameters();
//
//  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
//    {
//      truth_outputs[n] = 0.;
//      for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
//        truth_outputs[n] += get_rb_theta_expansion().eval_output_theta(n, q_l, mu)*
//          get_output_vector(n,q_l)->dot(*solution);
//    }
//
//
//      ExodusII_IO(get_mesh()).write_equation_systems ("truth.e",
//                                                      this->get_equation_systems());
//
//  // Get the X norm of the truth solution
//  // Useful for normalizing our true error data
//  inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//  Number truth_X_norm = std::sqrt(inner_product_storage_vector->dot(*solution));
//
//  return libmesh_real(truth_X_norm);
//
//  }


  unsigned int u_var;

};

///------------------------DWARFELEPHANTRBEVALUATION------------------------
class DwarfElephantRBEvaluationSteadyState : public RBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluationSteadyState(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem):
    RBEvaluation(comm),
//    offline_error_bound(false),
    fe_problem(fe_problem)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }

//    virtual Real get_stability_lower_bound()
//  {
//    _rb_scm_eval->set_parameters(get_parameters());
//    return _rb_scm_eval->get_SCM_LB();
//  }

  virtual Real get_stability_lower_bound()
  {
    const RBParameters & mu = get_parameters();

    Real min_mu = mu.get_value("mu_0");

    for (unsigned int  i = 1; i != mu.n_parameters(); i++)
    {
      const std::string mu_name = "mu_" + std::to_string(i);
      Real min_mu_i = std::min(min_mu, mu.get_value(mu_name));

      if (min_mu_i < min_mu)
        min_mu = min_mu_i;
    }

    return min_mu;
  }

  Real rb_solve(unsigned int N)
{
  LOG_SCOPE("rb_solve()", "RBEvaluation");

  if(N > get_n_basis_functions())
    libmesh_error_msg("ERROR: N cannot be larger than the number of basis functions in rb_solve");

  const RBParameters & mu = get_parameters();

  // Resize (and clear) the solution vector
  RB_solution.resize(N);

  // Assemble the RB system
  DenseMatrix<Number> RB_system_matrix(N,N);
  RB_system_matrix.zero();

  DenseMatrix<Number> RB_Aq_a;
  for(unsigned int q_a=0; q_a<_rb_theta_expansion.get_n_A_terms(); q_a++)
    {
      RB_Aq_vector[q_a].get_principal_submatrix(N, RB_Aq_a);

      RB_system_matrix.add(_rb_theta_expansion.eval_A_theta(q_a, mu), RB_Aq_a);
    }

  // Assemble the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  DenseVector<Number> RB_Fq_f;
  for(unsigned int q_f=0; q_f<_rb_theta_expansion.get_n_F_terms(); q_f++)
    {
      RB_Fq_vector[q_f].get_principal_subvector(N, RB_Fq_f);
      RB_rhs.add(_rb_theta_expansion.eval_F_theta(q_f, mu), RB_Fq_f);
    }

  // Solve the linear system
  if(N > 0)
    {
      RB_system_matrix.lu_solve(RB_rhs, RB_solution);
    }

  // Evaluate RB outputs
  DenseVector<Number> RB_output_vector_N;
  for(unsigned int n=0; n<_rb_theta_expansion.get_n_outputs(); n++)
    {
      RB_outputs[n] = 0.;
      for(unsigned int q_l=0; q_l<_rb_theta_expansion.get_n_output_terms(n); q_l++)
        {
          RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
          RB_outputs[n] += _rb_theta_expansion.eval_output_theta(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
        }
    }

  if(evaluate_RB_error_bound) // Calculate the error bounds
    {
      DwarfElephantRBConstructionSteadyState & sys_rb = fe_problem.es().get_system<DwarfElephantRBConstructionSteadyState>("RBSystem");
      // Evaluate the dual norm of the residual for RB_solution_vector

//      // slower but less error prone error bound (does not work in parallel)
//      epsilon_N = sys_rb.compute_residual_dual_norm(N);

      // faster but more error prone error bound (does work in parallel)
      epsilon_N = compute_residual_dual_norm(N);

      // Get lower bound for coercivity constant
      const Real alpha_LB = get_stability_lower_bound();
      // alpha_LB needs to be positive to get a valid error bound
      libmesh_assert_greater ( alpha_LB, 0. );

      // Evaluate the (absolute) error bound
      Real abs_error_bound = epsilon_N / residual_scaling_denom(alpha_LB);

      // Now compute the output error bounds
      for(unsigned int n=0; n<_rb_theta_expansion.get_n_outputs(); n++)
        {
          RB_output_error_bounds[n] = abs_error_bound * eval_output_dual_norm(n, mu);
        }

      return abs_error_bound;
    }
  else // Don't calculate the error bounds
    {
      // Just return -1. if we did not compute the error bound
      return -1.;
    }
}

  bool offline_error_bound;
  Real epsilon_N;
  FEProblemBase & fe_problem;
  RBP1T2EqualD2SteadyStateExpansion _rb_theta_expansion;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSESSTEADYSTATE_H
