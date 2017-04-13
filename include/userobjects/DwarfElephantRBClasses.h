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
#ifndef DWARFELEPHANTRBCLASSES_H
#define DWARFELEPHANTRBCLASSES_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/sparse_matrix.h"

// libMesh includes (RB package)
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"

///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "RBStructuresP1Theta3ThetaEqualMu.h"
#include "CacheBoundaries.h"

// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;

  class EquationSystems;
  class RBConstruction;
  class RBEvaluation;
}


///------------------------DWARFELEPHANTRBEVALUATION------------------------
class DwarfElephantRBEvaluation : public RBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluation(const libMesh::Parallel::Communicator & comm, CacheBoundaries * cache_boundaries):
    RBEvaluation(comm),
    cache_boundaries(cache_boundaries)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }

  virtual Real get_stability_lower_bound()
  {
//    const RBParameters & mu = get_parameters();
//    Real min_mu_1 = std::min(mu.get_value("mu_0"), mu.get_value("mu_1"));
//    return std::min(min_mu_1, mu.get_value("mu_2"));
    return 0.05;
  }

  CacheBoundaries & getBoundaryCache()
  {
    return *cache_boundaries;
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

//  libMesh::out << RB_Aq_vector[0] << std::endl;
  DenseMatrix<Number> RB_Aq_a;
  for(unsigned int q_a=0; q_a<_rb_theta_expansion.get_n_A_terms(); q_a++)
    {
      RB_Aq_vector[q_a].get_principal_submatrix(N, RB_Aq_a);

      RB_system_matrix.add(_rb_theta_expansion.eval_A_theta(q_a, mu), RB_Aq_a);
//      libMesh::out << RB_Aq_a << std::endl;
    }

//  libMesh::out << RB_system_matrix << std::endl;

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

Real compute_residual_dual_norm(const unsigned int N)
{
  LOG_SCOPE("compute_residual_dual_norm()", "RBEvaluation");

  const RBParameters & mu = get_parameters();

  // Use the stored representor inner product values
  // to evaluate the residual norm
  Number residual_norm_sq = 0.;

  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<_rb_theta_expansion.get_n_F_terms(); q_f1++)
    {
      for(unsigned int q_f2=q_f1; q_f2<_rb_theta_expansion.get_n_F_terms(); q_f2++)

        {
          Real delta = (q_f1==q_f2) ? 1. : 2.;
          residual_norm_sq += delta * libmesh_real(
                                                   _rb_theta_expansion.eval_F_theta(q_f1, mu)
                                                   * libmesh_conj(_rb_theta_expansion.eval_F_theta(q_f2, mu)) * Fq_representor_innerprods[q] );
//          libMesh::out << "F-F: " << residual_norm_sq << std::endl;
          q++;
        }
    }

  for(unsigned int q_f=0; q_f<_rb_theta_expansion.get_n_F_terms(); q_f++)
    {
      for(unsigned int q_a=0; q_a<_rb_theta_expansion.get_n_A_terms(); q_a++)
        {
          for(unsigned int i=0; i<N; i++)
            {
              Real delta = 2.;
              residual_norm_sq +=
                delta * libmesh_real( _rb_theta_expansion.eval_F_theta(q_f, mu) *
                                      libmesh_conj(_rb_theta_expansion.eval_A_theta(q_a, mu)) *
                                      libmesh_conj(RB_solution(i)) * Fq_Aq_representor_innerprods[q_f][q_a][i] );
//             libMesh::out << "F-A: " << residual_norm_sq << std::endl;
            }
        }
    }

  q=0;
  for(unsigned int q_a1=0; q_a1<_rb_theta_expansion.get_n_A_terms(); q_a1++)
    {
      for(unsigned int q_a2=q_a1; q_a2<_rb_theta_expansion.get_n_A_terms(); q_a2++)
        {
          Real delta = (q_a1==q_a2) ? 1. : 2.;

          for(unsigned int i=0; i<N; i++)
            {
              for(unsigned int j=0; j<N; j++)
                {
                  residual_norm_sq +=
                    delta * libmesh_real( libmesh_conj(_rb_theta_expansion.eval_A_theta(q_a1, mu)) *
                                          _rb_theta_expansion.eval_A_theta(q_a2, mu) *
                                          libmesh_conj(RB_solution(i)) * RB_solution(j) * Aq_Aq_representor_innerprods[q][i][j] );
//                libMesh::out << "A-A: " << residual_norm_sq << std::endl;
                }
            }

          q++;
        }
    }

  if(libmesh_real(residual_norm_sq) < 0.)
    {
      //    libMesh::out << "Warning: Square of residual norm is negative "
      //                 << "in RBSystem::compute_residual_dual_norm()" << std::endl;

      //     Sometimes this is negative due to rounding error,
      //     but when this occurs the error is on the order of 1.e-10,
      //     so shouldn't affect error bound much...
      residual_norm_sq = std::abs(residual_norm_sq);
    }

//  libMesh::out << "Final: " << residual_norm_sq << std::endl;

  return std::sqrt( libmesh_real(residual_norm_sq) );
}

  CacheBoundaries * cache_boundaries;
  RBP1Theta3ThetaEqualMuExpansion  _rb_theta_expansion;
};

///-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
class DwarfElephantRBConstruction : public RBConstruction
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBConstruction (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : Parent(es, name_in, number_in)
  {}

  // Destructor
  virtual ~DwarfElephantRBConstruction () { }

  // Type of the system
  typedef DwarfElephantRBConstruction _sys_type;

  // Type of the parent
  typedef RBConstruction Parent;

  // Initialize data structure
  virtual void init_data()
  {
    u_var = this->add_variable ("u", libMesh::FIRST);

    Parent::init_data();
  }

  void truth_assembly()
{
  LOG_SCOPE("truth_assembly()", "RBConstruction");

  const RBParameters & mu = get_parameters();

  RBEvaluation * _rb_eval = &get_rb_evaluation();
  DwarfElephantRBEvaluation * _new_rb_eval = dynamic_cast<DwarfElephantRBEvaluation *>(_rb_eval);
  CacheBoundaries & cache_boundaries = _new_rb_eval->getBoundaryCache();

  this->matrix->zero();
  this->rhs->zero();

  this->matrix->close();
  this->rhs->close();

  {
    // We should have already assembled the matrices
    // and vectors in the affine expansion, so
    // just use them

    for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
      {
        matrix->add(get_rb_theta_expansion().eval_A_theta(q_a, mu), *get_Aq(q_a));
//        get_Aq(q_a)->close();
//        cache_boundaries.resetBoundariesSubdomainStiffnessMatrix(*matrix,q_a);
//        matrix->close();
//        libMesh::out << *matrix << std::endl;
      }

    UniquePtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build(this->comm());
    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
      {
        *temp_vec = *get_Fq(q_f);
        temp_vec->scale( get_rb_theta_expansion().eval_F_theta(q_f, mu) );
        rhs->add(*temp_vec);
      }
  }

  this->matrix->close();
  this->rhs->close();
}

  Real truth_solve(int plot_solution)
{
  LOG_SCOPE("truth_solve()", "RBConstruction");

  truth_assembly();

  // truth_assembly assembles into matrix and rhs, so use those for the solve
  if (extra_linear_solver)
    {
      // If extra_linear_solver has been initialized, then we use it for the
      // truth solves.
      solve_for_matrix_and_rhs(*extra_linear_solver, *matrix, *rhs);

      if (assert_convergence)
        check_convergence(*extra_linear_solver);
    }
  else
    {
      solve_for_matrix_and_rhs(*get_linear_solver(), *matrix, *rhs);

      if (assert_convergence)
        check_convergence(*get_linear_solver());
    }



  const RBParameters & mu = get_parameters();

  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    {
      truth_outputs[n] = 0.;
      for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
        truth_outputs[n] += get_rb_theta_expansion().eval_output_theta(n, q_l, mu)*
          get_output_vector(n,q_l)->dot(*solution);
    }

     ExodusII_IO(get_mesh()).write_equation_systems ("truth.e",
                                                    this->get_equation_systems());


  // Get the X norm of the truth solution
  // Useful for normalizing our true error data
  inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
  Number truth_X_norm = std::sqrt(inner_product_storage_vector->dot(*solution));

//  libMesh::out << "Truth Norm: " << libmesh_real(truth_X_norm) << std::endl;

  return libmesh_real(truth_X_norm);
}

// overwrite the function to avoid negative values underneath the square root
  void enrich_RB_space()
  {
    LOG_SCOPE("enrich_RB_space()", "RBConstruction");
//
//    NumericVector<Number> * new_bf = NumericVector<Number>::build(this->comm()).release();
//    new_bf->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//    new_bf->zero();
//    new_bf->close();
//    libMesh::out << *new_bf << std::endl;


    NumericVector<Number> * new_bf = NumericVector<Number>::build(this->comm()).release();
    new_bf->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    *new_bf = *solution;

    for(unsigned int index=0; index<get_rb_evaluation().get_n_basis_functions(); index++)
    {
      inner_product_matrix->vector_mult(*inner_product_storage_vector, *new_bf);

      Number scalar =
        inner_product_storage_vector->dot(get_rb_evaluation().get_basis_function(index));
      new_bf->add(-scalar, get_rb_evaluation().get_basis_function(index));
//      libMesh::out << *new_bf << std::endl;
    }

    // Normalize new_bf
    inner_product_matrix->vector_mult(*inner_product_storage_vector, *new_bf);
    Number new_bf_norm = std::sqrt(std::abs(inner_product_storage_vector->dot(*new_bf)));

    if(new_bf_norm == 0.)
    {
      new_bf->zero(); // avoid potential nan's
    }
    else
    {
      new_bf->scale(1./new_bf_norm);
    }

    // load the new basis function into the basis_functions vector.
    get_rb_evaluation().basis_functions.push_back( new_bf );
}

void update_residual_terms(bool compute_inner_products = true)
{
  LOG_SCOPE("update_residual_terms()", "RBConstruction");

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
    {
      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
        {
          // Initialize the vector in which we'll store the representor
          if(!get_rb_evaluation().Aq_representor[q_a][i])
            {
              get_rb_evaluation().Aq_representor[q_a][i] = (NumericVector<Number>::build(this->comm()).release());
              get_rb_evaluation().Aq_representor[q_a][i]->init(this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
            }

          libmesh_assert(get_rb_evaluation().Aq_representor[q_a][i]->size()       == this->n_dofs()       &&
                         get_rb_evaluation().Aq_representor[q_a][i]->local_size() == this->n_local_dofs() );

          rhs->zero();
          get_Aq(q_a)->vector_mult(*rhs, get_rb_evaluation().get_basis_function(i));
          rhs->scale(-1.);
//          libMesh::out << get_Aq(0) << std::endl;

          solve_for_matrix_and_rhs(*inner_product_solver, *inner_product_matrix, *rhs);

          if (assert_convergence)
            check_convergence(*inner_product_solver);

          // Store the representor
          *get_rb_evaluation().Aq_representor[q_a][i] = *solution;
        }
    }

  // Now compute and store the inner products (if requested)
  if (compute_inner_products)
    {

      for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
        {
          inner_product_matrix->vector_mult(*inner_product_storage_vector,*Fq_representor[q_f]);

          for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
            {
              for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
                {
                  get_rb_evaluation().Fq_Aq_representor_innerprods[q_f][q_a][i] =
                    inner_product_storage_vector->dot(*get_rb_evaluation().Aq_representor[q_a][i]);
                }
            }
        }

      unsigned int q=0;
      for(unsigned int q_a1=0; q_a1<get_rb_theta_expansion().get_n_A_terms(); q_a1++)
        {
          for(unsigned int q_a2=q_a1; q_a2<get_rb_theta_expansion().get_n_A_terms(); q_a2++)
            {
              for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
                {
                  for(unsigned int j=0; j<RB_size; j++)
                    {
                      inner_product_matrix->vector_mult(*inner_product_storage_vector, *get_rb_evaluation().Aq_representor[q_a2][j]);
                      get_rb_evaluation().Aq_Aq_representor_innerprods[q][i][j] =
                        inner_product_storage_vector->dot(*get_rb_evaluation().Aq_representor[q_a1][i]);

                      if(i != j)
                        {
                          inner_product_matrix->vector_mult(*inner_product_storage_vector, *get_rb_evaluation().Aq_representor[q_a2][i]);
                          get_rb_evaluation().Aq_Aq_representor_innerprods[q][j][i] =
                            inner_product_storage_vector->dot(*get_rb_evaluation().Aq_representor[q_a1][j]);
                        }
                    }
                }
              q++;
            }
        }
    } // end if (compute_inner_products)
}

  unsigned int u_var;
  DwarfElephantRBEvaluation * _new_rb_eval;

  friend class DwarfElephantOfflineStage;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSES_H
