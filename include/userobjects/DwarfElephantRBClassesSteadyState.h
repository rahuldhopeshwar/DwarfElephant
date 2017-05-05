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
 * NOTE: ENSURE THAT THE CLASS IS INHERITING FROM THE CORRECT RBSTRUCTURES.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBCLASSESSTEADYSTATE_H
#define DWARFELEPHANTRBCLASSESSTEADYSTATE_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"

// libMesh includes (RB package)
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"

///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "CacheBoundaries.h"
#include "RBStructuresP1Theta3ThetaEqualMuSteadyState.h"
#include "RBStructuresP1Theta5ThetaEqualMuSteadyState.h"

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

///------------------------DWARFELEPHANTRBEVALUATION------------------------
class DwarfElephantRBEvaluationSteadyState : public RBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluationSteadyState(const libMesh::Parallel::Communicator & comm): //, CacheBoundaries * cache_boundaries, FEProblemBase & fe_problem):
    RBEvaluation(comm)
//    cache_boundaries(cache_boundaries),
//    fe_problem(fe_problem)
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

//  CacheBoundaries & getBoundaryCache()
//  {
//    return *cache_boundaries;
//  }
//
//  FEProblemBase & getFEProblemBase()
//  {
//    return fe_problem;
//  }
//
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
      // Evaluate the dual norm of the residual for RB_solution_vector
      Real epsilon_N = compute_residual_dual_norm(N);
//      libMesh::out << "epsilon_N: " << epsilon_N << std::endl;

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
//
//Real compute_residual_dual_norm(const unsigned int N)
//{
//  LOG_SCOPE("compute_residual_dual_norm()", "RBEvaluation");
//
//  const RBParameters & mu = get_parameters();
//
//  // Use the stored representor inner product values
//  // to evaluate the residual norm
//  Number residual_norm_sq = 0.;
//
//  unsigned int q=0;
//  for(unsigned int q_f1=0; q_f1<_rb_theta_expansion.get_n_F_terms(); q_f1++)
//    {
//      for(unsigned int q_f2=q_f1; q_f2<_rb_theta_expansion.get_n_F_terms(); q_f2++)
//        {
//          Real delta = (q_f1==q_f2) ? 1. : 2.;
//          residual_norm_sq += delta * libmesh_real(
//                                                   _rb_theta_expansion.eval_F_theta(q_f1, mu)
//                                                   * libmesh_conj(_rb_theta_expansion.eval_F_theta(q_f2, mu)) * Fq_representor_innerprods[q] );
////          libMesh::out << "Fq_representor_inner: " << Fq_representor_innerprods[q] << std::endl;
//          q++;
////          libMesh::out << "Norm 1:" << residual_norm_sq << std::endl;
//        }
//    }
//
//  for(unsigned int q_f=0; q_f<_rb_theta_expansion.get_n_F_terms(); q_f++)
//    {
//      for(unsigned int q_a=0; q_a<_rb_theta_expansion.get_n_A_terms(); q_a++)
//        {
//          for(unsigned int i=0; i<N; i++)
//            {
//              Real delta = 2.;
//              residual_norm_sq +=
//                delta * libmesh_real( _rb_theta_expansion.eval_F_theta(q_f, mu) *
//                                      libmesh_conj(_rb_theta_expansion.eval_A_theta(q_a, mu)) *
//                                      libmesh_conj(RB_solution(i)) * Fq_Aq_representor_innerprods[q_f][q_a][i] );
////                libMesh::out << "Norm 2:" << residual_norm_sq << std::endl;
//            }
//        }
//    }
//
//  q=0;
//  for(unsigned int q_a1=0; q_a1<_rb_theta_expansion.get_n_A_terms(); q_a1++)
//    {
//      for(unsigned int q_a2=q_a1; q_a2<_rb_theta_expansion.get_n_A_terms(); q_a2++)
//        {
//          Real delta = (q_a1==q_a2) ? 1. : 2.;
//
//          for(unsigned int i=0; i<N; i++)
//            {
//              for(unsigned int j=0; j<N; j++)
//                {
//                  residual_norm_sq +=
//                    delta * libmesh_real( libmesh_conj(_rb_theta_expansion.eval_A_theta(q_a1, mu)) *
//                                          _rb_theta_expansion.eval_A_theta(q_a2, mu) *
//                                          libmesh_conj(RB_solution(i)) * RB_solution(j) * Aq_Aq_representor_innerprods[q][i][j] );
//                }
//            }
//
//          q++;
////          libMesh::out << "Norm 3:" << residual_norm_sq << std::endl;
//        }
//    }
//
//  if(libmesh_real(residual_norm_sq) < 0.)
//    {
//      //    libMesh::out << "Warning: Square of residual norm is negative "
//      //                 << "in RBSystem::compute_residual_dual_norm()" << std::endl;
//
//      //     Sometimes this is negative due to rounding error,
//      //     but when this occurs the error is on the order of 1.e-10,
//      //     so shouldn't affect error bound much...
//      residual_norm_sq = std::abs(residual_norm_sq);
//    }
//
////  libMesh::out << "Norm total:" << residual_norm_sq << std::endl;
////  libMesh::out << "Norm sqrt:" << std::sqrt( libmesh_real(residual_norm_sq)) << std::endl;
//  return std::sqrt( libmesh_real(residual_norm_sq) );
//}

//  CacheBoundaries * cache_boundaries;
//  FEProblemBase & fe_problem;
  RBP1Theta3ThetaEqualMuExpansionSteadyState _rb_theta_expansion;
};

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

//  Real get_RB_error_bound()
//{
//  get_rb_evaluation().set_parameters( get_parameters() );
//
//  Real error_bound = get_rb_evaluation().rb_solve(get_rb_evaluation().get_n_basis_functions());
//
////  if (normalize_rb_bound_in_greedy)
////    {
////      Real error_bound_normalization = get_rb_evaluation().get_error_bound_normalization();
////
////      if ((error_bound < abs_training_tolerance) ||
////          (error_bound_normalization < abs_training_tolerance))
////        {
////          // We don't want to normalize this error bound if the bound or the
////          // normalization value are below the absolute tolerance. Hence do nothing
////          // in this case.
////        }
////      else
////        error_bound /= error_bound_normalization;
////    }
//
//  libMesh::out << "Internal error bound: " << error_bound << std::endl;
//  return error_bound;
//}

//  Real train_reduced_basis(const bool resize_rb_eval_data=true)
//{
//  LOG_SCOPE("train_reduced_basis()", "RBConstruction");
//
//  int count = 0;
//
//  // initialize rb_eval's parameters
//  get_rb_evaluation().initialize_parameters(*this);
//
//  // possibly resize data structures according to Nmax
//  if(resize_rb_eval_data)
//    {
//      get_rb_evaluation().resize_data_structures(get_Nmax());
//    }
//
//  // Clear the Greedy param list
//  for (std::size_t i=0; i<get_rb_evaluation().greedy_param_list.size(); i++)
//    get_rb_evaluation().greedy_param_list[i].clear();
//
//  get_rb_evaluation().greedy_param_list.clear();
//
//  Real training_greedy_error;
//
//
//  // If we are continuing from a previous training run,
//  // we might already be at the max number of basis functions.
//  // If so, we can just return.
//  if(get_rb_evaluation().get_n_basis_functions() >= get_Nmax())
//    {
//      libMesh::out << "Maximum number of basis functions reached: Nmax = "
//                   << get_Nmax() << std::endl;
//      return 0.;
//    }
//
//
//  // Compute the dual norms of the outputs if we haven't already done so
//  compute_output_dual_innerprods();
//
//  // Compute the Fq Riesz representor dual norms if we haven't already done so
//  compute_Fq_representor_innerprods();
//
//  libMesh::out << std::endl << "---- Performing Greedy basis enrichment ----" << std::endl;
//  Real initial_greedy_error = 0.;
//  bool initial_greedy_error_initialized = false;
//  while(true)
//    {
//      libMesh::out << std::endl << "---- Basis dimension: "
//                   << get_rb_evaluation().get_n_basis_functions() << " ----" << std::endl;
//
//      if( count > 0 || (count==0 && use_empty_rb_solve_in_greedy) )
//        {
//          libMesh::out << "Performing RB solves on training set" << std::endl;
//          training_greedy_error = compute_max_error_bound();
//
//          libMesh::out << "Maximum error bound is " << training_greedy_error << std::endl << std::endl;
//
//          // record the initial error
//          if (!initial_greedy_error_initialized)
//            {
//              initial_greedy_error = training_greedy_error;
//              initial_greedy_error_initialized = true;
//            }
//
//          // Break out of training phase if we have reached Nmax
//          // or if the training_tolerance is satisfied.
//          if (greedy_termination_test(training_greedy_error, initial_greedy_error, count))
//            break;
//        }
//
//      libMesh::out << "Performing truth solve at parameter:" << std::endl;
//      print_parameters();
//
//      // Update the list of Greedily selected parameters
//      this->update_greedy_param_list();
//
//      // Perform an Offline truth solve for the current parameter
//      truth_solve(-1);
//
//      // Add orthogonal part of the snapshot to the RB space
//      libMesh::out << "Enriching the RB space" << std::endl;
//      enrich_RB_space();
//
//      update_system();
//
//      // Increment counter
//      count++;
//    }
//  this->update_greedy_param_list();
//
//  return training_greedy_error;
//}
//
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
//    ExodusII_IO(get_mesh()).write_equation_systems ("truth.e",
//                                                      this->get_equation_systems());
//  if(plot_solution > 0)
//    {
//#if defined(LIBMESH_USE_COMPLEX_NUMBERS)
//      GMVIO(get_mesh()).write_equation_systems ("truth.gmv",
//                                                this->get_equation_systems());
//#else
//#ifdef LIBMESH_HAVE_EXODUS_API
//      ExodusII_IO(get_mesh()).write_equation_systems ("truth.e",
//                                                      this->get_equation_systems());
//#endif
//#endif
//    }
//
//  // Get the X norm of the truth solution
//  // Useful for normalizing our true error data
//  inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//  Number truth_X_norm = std::sqrt(inner_product_storage_vector->dot(*solution));
//
//  return libmesh_real(truth_X_norm);
//}
//
//void truth_assembly()
//{
//  LOG_SCOPE("truth_assembly()", "RBConstruction");
//
//  const RBParameters & mu = get_parameters();
//
//  this->matrix->zero();
//  this->rhs->zero();
//
//  this->matrix->close();
//  this->rhs->close();
//
//  {
//    // We should have already assembled the matrices
//    // and vectors in the affine expansion, so
//    // just use them
//
//   RBEvaluation * rb_eval = &get_rb_evaluation();
//   DwarfElephantRBEvaluation * rb_eval_dwarf = dynamic_cast<DwarfElephantRBEvaluation *>(rb_eval);
//
//    for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
//      {
//        matrix->add(get_rb_theta_expansion().eval_A_theta(q_a, mu), *get_Aq(q_a));
//      }
//
//    UniquePtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build(this->comm());
//    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//    for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
//      {
//        *temp_vec = *get_Fq(q_f);
//        temp_vec->scale( get_rb_theta_expansion().eval_F_theta(q_f, mu) );
//        rhs->add(*temp_vec);
//      }
//  }
//
//  this->matrix->close();
//  this->rhs->close();
//}
//
//void compute_output_dual_innerprods()
//{
//  // Skip calculations if we've already computed the output dual norms
//  if(!output_dual_innerprods_computed)
//    {
//      // Short circuit if we don't have any outputs
//      if( get_rb_theta_expansion().get_n_outputs() == 0 )
//        {
//          output_dual_innerprods_computed = true;
//          return;
//        }
//
//      // Only log if we get to here
//      LOG_SCOPE("compute_output_dual_innerprods()", "RBConstruction");
//
//      libMesh::out << "Compute output dual inner products" << std::endl;
//
//      // Find out the largest value of Q_l
//      unsigned int max_Q_l = 0;
//      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
//        max_Q_l = (get_rb_theta_expansion().get_n_output_terms(n) > max_Q_l) ? get_rb_theta_expansion().get_n_output_terms(n) : max_Q_l;
//
//      std::vector< NumericVector<Number> * > L_q_representor(max_Q_l);
//      for(unsigned int q=0; q<max_Q_l; q++)
//        {
//          L_q_representor[q] = (NumericVector<Number>::build(this->comm()).release());
//          L_q_representor[q]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//        }
//
//      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
//        {
//          for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
//            {
//              rhs->zero();
//              rhs->add(1., *get_output_vector(n,q_l));
//
//              // Use the main linear solver here instead of the inner_product solver, since
//              // get_matrix_for_output_dual_solves() may not return the inner product matrix.
//              solve_for_matrix_and_rhs(*get_linear_solver(), get_matrix_for_output_dual_solves(), *rhs);
//
//              // We possibly perform multiple solves here with the same matrix, hence
//              // set reuse_preconditioner(true) (and set it back to false again below
//              // at the end of this function).
//              linear_solver->reuse_preconditioner(true);
//
//              if (assert_convergence)
//                check_convergence(*get_linear_solver());
//
//              *L_q_representor[q_l] = *solution;
//            }
//
//          unsigned int q=0;
//          for(unsigned int q_l1=0; q_l1<get_rb_theta_expansion().get_n_output_terms(n); q_l1++)
//            {
//              get_matrix_for_output_dual_solves().vector_mult(*inner_product_storage_vector, *L_q_representor[q_l1]);
//
//              for(unsigned int q_l2=q_l1; q_l2<get_rb_theta_expansion().get_n_output_terms(n); q_l2++)
//                {
//                  output_dual_innerprods[n][q] = L_q_representor[q_l2]->dot(*inner_product_storage_vector);
//                  libMesh::out << "output_dual_innerprods[" << n << "][" << q << "] = " << output_dual_innerprods[n][q] << std::endl;
//
//                  q++;
//                }
//            }
//        }
//
//      // Finally clear the L_q_representor vectors
//      for(unsigned int q=0; q<max_Q_l; q++)
//        {
//          if(L_q_representor[q])
//            {
//              delete L_q_representor[q];
//              L_q_representor[q] = libmesh_nullptr;
//            }
//        }
//
//      // We may not need to use linear_solver again (e.g. this would happen if we use
//      // extra_linear_solver for the truth_solves). As a result, let's clear linear_solver
//      // to release any memory it may be taking up. If we do need it again, it will
//      // be initialized when necessary.
//      linear_solver->clear();
//      linear_solver->reuse_preconditioner(false);
//
//      output_dual_innerprods_computed = true;
//    }
//
//  get_rb_evaluation().output_dual_innerprods = output_dual_innerprods;
//}
//
//// overwrite the function to avoid negative values underneath the square root
//  void enrich_RB_space()
//  {
//    LOG_SCOPE("enrich_RB_space()", "RBConstruction");
//
//    NumericVector<Number> * new_bf = NumericVector<Number>::build(this->comm()).release();
//    new_bf->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
////    solution->scale(-1);
//    *new_bf = *solution;
//
//    for(unsigned int index=0; index<get_rb_evaluation().get_n_basis_functions(); index++)
//    {
//      inner_product_matrix->vector_mult(*inner_product_storage_vector, *new_bf);
//
//      Number scalar =
//        inner_product_storage_vector->dot(get_rb_evaluation().get_basis_function(index));
//      new_bf->add(-scalar, get_rb_evaluation().get_basis_function(index));
//    }
//
//    // Normalize new_bf
//    inner_product_matrix->vector_mult(*inner_product_storage_vector, *new_bf);
//    Number new_bf_norm = std::sqrt(std::abs(inner_product_storage_vector->dot(*new_bf)) );
//
//    libMesh::out << "RB Norm: " << new_bf_norm << std::endl;
//
//    if(new_bf_norm == 0.)
//    {
//      new_bf->zero(); // avoid potential nan's
//    }
//    else
//    {
//      new_bf->scale(1./new_bf_norm);
//    }
//
//    // load the new basis function into the basis_functions vector.
//    get_rb_evaluation().basis_functions.push_back( new_bf );
//}
//
//Real compute_max_error_bound()
//{
//  LOG_SCOPE("compute_max_error_bound()", "RBConstruction");
//
//  // Treat the case with no parameters in a special way
//  if(get_n_params() == 0)
//    {
//      Real max_val;
//      if(std::numeric_limits<Real>::has_infinity)
//        {
//          max_val = std::numeric_limits<Real>::infinity();
//        }
//      else
//        {
//          max_val = std::numeric_limits<Real>::max();
//        }
//
//      // Make sure we do at least one solve, but otherwise return a zero error bound
//      // when we have no parameters
//      return (get_rb_evaluation().get_n_basis_functions() == 0) ? max_val : 0.;
//    }
//
//  training_error_bounds.resize(this->get_local_n_training_samples());
//
//  // keep track of the maximum error
//  unsigned int max_err_index = 0;
//  Real max_err = 0.;
//
//  numeric_index_type first_index = get_first_local_training_index();
//  for(unsigned int i=0; i<get_local_n_training_samples(); i++)
//    {
//      // Load training parameter i, this is only loaded
//      // locally since the RB solves are local.
//      set_params_from_training_set( first_index+i );
//
//      training_error_bounds[i] = get_RB_error_bound();
////      libMesh::out << "Training error: " << training_error_bounds[i] << std::endl;
//
//      if(training_error_bounds[i] > max_err)
//        {
//          max_err_index = i;
//          max_err = training_error_bounds[i];
//        }
//    }
//
//  std::pair<numeric_index_type, Real> error_pair(first_index+max_err_index, max_err);
//  get_global_max_error_pair(this->comm(),error_pair);
//
//  // If we have a serial training set (i.e. a training set that is the same on all processors)
//  // just set the parameters on all processors
//  if(serial_training_set)
//    {
//      set_params_from_training_set( error_pair.first );
//    }
//  // otherwise, broadcast the parameter that produced the maximum error
//  else
//    {
//      unsigned int root_id=0;
//      if( (get_first_local_training_index() <= error_pair.first) &&
//          (error_pair.first < get_last_local_training_index()) )
//        {
//          set_params_from_training_set( error_pair.first );
//          root_id = this->processor_id();
//        }
//
//      this->comm().sum(root_id); // root_id is only non-zero on one processor
//      broadcast_parameters(root_id);
//    }
//
//  return error_pair.second;
//}
//
//Real get_RB_error_bound()
//{
//  get_rb_evaluation().set_parameters( get_parameters() );
//
//  Real error_bound = get_rb_evaluation().rb_solve(get_rb_evaluation().get_n_basis_functions());
////  libMesh::out << "RB error bound: " << error_bound << std::endl;
//
//  return error_bound;
//}
//
//void compute_Fq_representor_innerprods(bool compute_inner_products=true)
//{
//
//  // Skip calculations if we've already computed the Fq_representors
//  if(!Fq_representor_innerprods_computed)
//    {
//      // Only log if we get to here
//      LOG_SCOPE("compute_Fq_representor_innerprods()", "RBConstruction");
//
//                RBEvaluation * rb_eval = &get_rb_evaluation();
//          DwarfElephantRBEvaluation * rb_eval_dwarf = dynamic_cast<DwarfElephantRBEvaluation *>(rb_eval);
//          NumericVector<Number> * Fq_representor_total = (NumericVector<Number>::build(this->comm()).release());
//          Fq_representor_total->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//          Fq_representor_total->zero();
//
//      for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
//        {
//          if(!Fq_representor[q_f])
//            {
//              Fq_representor[q_f] = (NumericVector<Number>::build(this->comm()).release());
//              Fq_representor[q_f]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//            }
//
//          libmesh_assert(Fq_representor[q_f]->size()       == this->n_dofs()       &&
//                         Fq_representor[q_f]->local_size() == this->n_local_dofs() );
//
//          rhs->zero();
//          rhs->add(1., *get_Fq(q_f));
//
//
//          solve_for_matrix_and_rhs(*inner_product_solver, *inner_product_matrix, *rhs);
//
//          if (assert_convergence)
//            check_convergence(*inner_product_solver);
//
//          *Fq_representor[q_f] = *solution;
//
//        }
//
//      if (compute_inner_products)
//        {
//          unsigned int q=0;
//
//          for(unsigned int q_f1=0; q_f1<get_rb_theta_expansion().get_n_F_terms(); q_f1++)
//            {
//            inner_product_matrix->vector_mult(*inner_product_storage_vector, *Fq_representor[q_f1]);
//
//              for(unsigned int q_f2=q_f1; q_f2<get_rb_theta_expansion().get_n_F_terms(); q_f2++)
//                {
//                  Fq_representor_innerprods[q] = inner_product_storage_vector->dot(*Fq_representor[q_f2]);
//                  q++;
//                }
//            }
//        } // end if (compute_inner_products)
//
//      Fq_representor_innerprods_computed = true;
//    }
//
//  get_rb_evaluation().Fq_representor_innerprods = Fq_representor_innerprods;
//}



  unsigned int u_var;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSESSTEADYSTATE_H
