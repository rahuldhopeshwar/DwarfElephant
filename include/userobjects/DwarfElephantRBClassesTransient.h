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
#ifndef DWARFELEPHANTRBCLASSESTRANSIENT_H
#define DWARFELEPHANTRBCLASSESTRANSIENT_H

///---------------------------------INCLUDES--------------------------------
//#if defined(LIBMESH_HAVE_SLEPC) && defined(LIBMESH_HAVE_GLPK)

// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"

// libMesh includes (RB package)
#include "libmesh/transient_rb_evaluation.h"
#include "libmesh/transient_rb_construction.h"

///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBStructuresP1T1EqualF1O1Transient.h"
#include "DwarfElephantRBStructuresP1T2EqualF1O1Transient.h"
#include "DwarfElephantRBStructuresP1T3EqualF1O1Transient.h"
#include "DwarfElephantRBStructuresP1T3EqualF3O1Transient.h"
#include "DwarfElephantRBStructuresP1T4EqualF1O1Transient.h"
#include "DwarfElephantRBStructuresP1T5EqualF1O1Transient.h"
#include "DwarfElephantRBStructuresP1T5EqualF3O1Transient.h"

#include "FEProblemBase.h"

// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;

  class EquationSystems;
  class TransientRBConstruction;
  class TransientRBEvaluation;
}

//using libMesh::RBSCMConstruction;
//using libMesh::RBSCMEvaluation;

///------------------------DWARFELEPHANTRBEVALUATION------------------------
class DwarfElephantRBEvaluationTransient : public TransientRBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluationTransient(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem):
    TransientRBEvaluation(comm),
    fe_problem(fe_problem)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }

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
  LOG_SCOPE("rb_solve()", "TransientRBEvaluation");

  if (N > get_n_basis_functions())
    libmesh_error_msg("ERROR: N cannot be larger than the number of basis functions in rb_solve");

  const RBParameters & mu = get_parameters();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();
  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  const unsigned int n_time_steps = get_n_time_steps();
  const Real dt                   = get_delta_t();
  const Real euler_theta          = get_euler_theta();

  // Resize the RB and error bound vectors
  error_bound_all_k.resize(n_time_steps+1);
  RB_outputs_all_k.resize(trans_theta_expansion.get_n_outputs());
  RB_output_error_bounds_all_k.resize(trans_theta_expansion.get_n_outputs());
  for (unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
    {
      RB_outputs_all_k[n].resize(n_time_steps+1, 0.);
      RB_output_error_bounds_all_k[n].resize(n_time_steps+1, 0.);
    }

  // First assemble the mass matrix
  DenseMatrix<Number> RB_mass_matrix_N(N,N);
  RB_mass_matrix_N.zero();
  DenseMatrix<Number> RB_M_q_m;
  for (unsigned int q_m=0; q_m<Q_m; q_m++)
    {
      RB_M_q_vector[q_m].get_principal_submatrix(N, RB_M_q_m);
      RB_mass_matrix_N.add(trans_theta_expansion.eval_M_theta(q_m, mu), RB_M_q_m);
    }

  RB_LHS_matrix.resize(N,N);
  RB_LHS_matrix.zero();

  RB_RHS_matrix.resize(N,N);
  RB_RHS_matrix.zero();

  RB_LHS_matrix.add(1./dt, RB_mass_matrix_N);
  RB_RHS_matrix.add(1./dt, RB_mass_matrix_N);

  DenseMatrix<Number> RB_Aq_a;
  for (unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      RB_Aq_vector[q_a].get_principal_submatrix(N, RB_Aq_a);

      RB_LHS_matrix.add(       euler_theta*trans_theta_expansion.eval_A_theta(q_a,mu), RB_Aq_a);
      RB_RHS_matrix.add( -(1.-euler_theta)*trans_theta_expansion.eval_A_theta(q_a,mu), RB_Aq_a);
    }

  // Add forcing terms
  DenseVector<Number> RB_Fq_f;
  RB_RHS_save.resize(N);
  RB_RHS_save.zero();
  for (unsigned int q_f=0; q_f<Q_f; q_f++)
    {
      RB_Fq_vector[q_f].get_principal_subvector(N, RB_Fq_f);
      RB_RHS_save.add(trans_theta_expansion.eval_F_theta(q_f,mu), RB_Fq_f);
    }

  // Set system time level to 0
  set_time_step(0);

  // Resize/clear the solution vector
  RB_solution.resize(N);

  // Load the initial condition into RB_solution
  if (N > 0)
    {
      RB_solution = RB_initial_condition_all_N[N-1];
    }

  // Resize/clear the old solution vector
  old_RB_solution.resize(N);

  // Initialize the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  // Initialize the vectors storing solution data
  RB_temporal_solution_data.resize(n_time_steps+1);
  for (unsigned int time_level=0; time_level<=n_time_steps; time_level++)
    {
      RB_temporal_solution_data[time_level].resize(N);
    }
  // and load the initial data
  RB_temporal_solution_data[0] = RB_solution;

  // Set outputs at initial time
  {
    DenseVector<Number> RB_output_vector_N;
    for (unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
      {
        RB_outputs_all_k[n][0] = 0.;
        for (unsigned int q_l=0; q_l<trans_theta_expansion.get_n_output_terms(n); q_l++)
          {
            RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
            RB_outputs_all_k[n][0] += trans_theta_expansion.eval_output_theta(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
          }
      }
  }

  // Initialize error bounds, if necessary
  Real error_bound_sum = 0.;
  Real alpha_LB = 0.;
  if (evaluate_RB_error_bound)
    {
      if (N > 0)
        {
          error_bound_sum += pow( initial_L2_error_all_N[N-1], 2.);
        }

      // Set error bound at the initial time
      error_bound_all_k[get_time_step()] = std::sqrt(error_bound_sum);

      // Compute the outputs and associated error bounds at the initial time
      DenseVector<Number> RB_output_vector_N;
      for (unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
        {
          RB_outputs_all_k[n][0] = 0.;
          for (unsigned int q_l=0; q_l<trans_theta_expansion.get_n_output_terms(n); q_l++)
            {
              RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
              RB_outputs_all_k[n][0] += trans_theta_expansion.eval_output_theta(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
            }

          RB_output_error_bounds_all_k[n][0] = error_bound_all_k[0] * eval_output_dual_norm(n,mu);
        }

      alpha_LB = get_stability_lower_bound();

      // Precompute time-invariant parts of the dual norm of the residual.
      cache_online_residual_terms(N);
    }

  for (unsigned int time_level=1; time_level<=n_time_steps; time_level++)
    {
      set_time_step(time_level);
      old_RB_solution = RB_solution;

      // Compute RB_rhs, as RB_LHS_matrix x old_RB_solution
      RB_RHS_matrix.vector_mult(RB_rhs, old_RB_solution);

      // Add forcing terms
      RB_rhs.add(get_control(time_level), RB_RHS_save);

      if (N > 0)
        {
          RB_LHS_matrix.lu_solve(RB_rhs, RB_solution);
        }

      // Save RB_solution for current time level
      RB_temporal_solution_data[time_level] = RB_solution;

      // Evaluate outputs
      DenseVector<Number> RB_output_vector_N;
      for (unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
        {
          RB_outputs_all_k[n][time_level] = 0.;
          for (unsigned int q_l=0; q_l<trans_theta_expansion.get_n_output_terms(n); q_l++)
            {
              RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
              RB_outputs_all_k[n][time_level] += trans_theta_expansion.eval_output_theta(n,q_l,mu)*
                RB_output_vector_N.dot(RB_solution);
            }
        }

      // Calculate RB error bounds
      if (evaluate_RB_error_bound)
        {
          // Evaluate the dual norm of the residual for RB_solution_vector
          // Real epsilon_N = uncached_compute_residual_dual_norm(N);
          Real epsilon_N = compute_residual_dual_norm(N);

          error_bound_sum += residual_scaling_numer(alpha_LB) * pow(epsilon_N, 2.);

          // store error bound at time-level _k
          error_bound_all_k[time_level] = std::sqrt(error_bound_sum/residual_scaling_denom(alpha_LB));

          // Now evaluated output error bounds
          for (unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
            {
              RB_output_error_bounds_all_k[n][time_level] = error_bound_all_k[time_level] *
                eval_output_dual_norm(n,mu);
            }
        }
    }

  _rb_solve_data_cached = true ;

  if (evaluate_RB_error_bound) // Calculate the error bounds
    {
      return error_bound_all_k[n_time_steps];
    }
  else // Don't calculate the error bounds
    {
      // Just return -1. if we did not compute the error bound
      return -1.;
    }
}

  FEProblemBase & get_fe_problem() {return fe_problem;}

  FEProblemBase & fe_problem;
  DwarfElephantRBP1T5EqualF3O1TransientExpansion _rb_theta_expansion;
};


///-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
class DwarfElephantRBConstructionTransient : public TransientRBConstruction
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBConstructionTransient (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : Parent(es, name_in, number_in)
  {}

  // Destructor
  virtual ~DwarfElephantRBConstructionTransient () { }

  // Type of the system
  typedef DwarfElephantRBConstructionTransient _sys_type;

  // Type of the parent
  typedef TransientRBConstruction Parent;

  // Initialize data structure
  virtual void init_data()
  {
    u_var = this->add_variable (get_equation_systems().get_system(0).variable_name(0) + "(RB)", libMesh::FIRST);

    Parent::init_data();
  }

//  void initialize_truth ()
//{
//  DwarfElephantRBEvaluationTransient & _trans_rb_eval = cast_ref<DwarfElephantRBEvaluationTransient &>(get_rb_evaluation());
//  *this->solution = *_trans_rb_eval.get_fe_problem().es().get_system<TransientNonlinearImplicitSystem>("rb0").solution;
//  this->solution->close();
//  this->update();
//}

//  Real compute_residual_dual_norm(const unsigned int N)
//{
//   LOG_SCOPE("compute_residual_dual_norm()", "RBConstruction");
//
//   // Put the residual in rhs in order to compute the norm of the Riesz representor
//   // Note that this only works in serial since otherwise each processor will
//   // have a different parameter value during the Greedy training.
//
//   UniquePtr< NumericVector<Number> > RB_sol = NumericVector<Number>::build(this->comm());
//   RB_sol->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//
//   UniquePtr< NumericVector<Number> > temp = NumericVector<Number>::build(this->comm());
//   temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//
//   for(unsigned int i=0; i<N; i++)
//   {
//     RB_sol->add(get_rb_evaluation().RB_solution(i), get_rb_evaluation().get_basis_function(i));
//   }
//
//   this->truth_assembly();
//   matrix->vector_mult(*temp, *RB_sol);
//   rhs->add(-1., *temp);
//
//   // Then solve to get the Reisz representor
//   matrix->zero();
//   matrix->add(1., *inner_product_matrix);
////   if(constrained_problem)
////     matrix->add(1., *constraint_matrix);
//
//   solution->zero();
//   solve();
//   // Make sure we didn't max out the number of iterations
//   if( (this->n_linear_iterations() >=
//        this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
//       (this->final_linear_residual() >
//        this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
//   {
//     libmesh_error_msg("Warning: Linear solver may not have converged! Final linear residual = "
//                       << this->final_linear_residual() << ", number of iterations = "
//                       << this->n_linear_iterations());
//   }
//
//   inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//
//   Real slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);
//
//   return std::sqrt( libmesh_real(slow_residual_norm_sq) );
//}

  unsigned int u_var;

};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSESTRANSIENT_H
