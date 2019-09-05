/**
 * In this class subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
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

//include additional part for constraint_matrix
//#include "libmesh/coupling_matrix.h"

// MOOSE includes
#include "FEProblemBase.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemSteadyState.h"

// #include "DwarfElephantRBStructuresT1F1O1SteadyState.h"
// #include "DwarfElephantRBStructuresT2F1O1SteadyState.h"
// #include "DwarfElephantRBStructuresT2F1O3SteadyState.h"
// #include "DwarfElephantRBStructuresT2F1O5SteadyState.h"
// #include "DwarfElephantRBStructuresT2F3O1SteadyState.h"
// #include "DwarfElephantRBStructuresT2F2O10SteadyState.h"
// #include "DwarfElephantRBStructuresT2F3O10SteadyState.h"
// #include "DwarfElephantRBStructuresT3F1O1SteadyState.h"
// #include "DwarfElephantRBStructuresT3F2O1SteadyState.h"
// #include "DwarfElephantRBStructuresT3F3O1SteadyState.h"
// #include "DwarfElephantRBStructuresT3F4O1SteadyState.h"
// #include "DwarfElephantRBStructuresT3F1O3SteadyState.h"
// #include "DwarfElephantRBStructuresT3F3O10SteadyState.h"
// #include "DwarfElephantRBStructuresT3F4O2347SteadyState.h"
// #include "DwarfElephantRBStructuresT4F1O1SteadyState.h"
// #include "DwarfElephantRBStructuresT4F2O1SteadyState.h"
// #include "DwarfElephantRBStructuresT4F4O1SteadyState.h"
// #include "DwarfElephantRBStructuresT4F7O1SteadyState.h"
// #include "DwarfElephantRBStructuresT4F10O1SteadyState.h"
// #include "DwarfElephantRBStructuresT4F5O10SteadyState.h"
// #include "DwarfElephantRBStructuresT4F1O32SteadyState.h"
// #include "DwarfElephantRBStructuresT5F1O1SteadyState.h"
// #include "DwarfElephantRBStructuresT5F3O1SteadyState.h"
// #include "DwarfElephantRBStructuresT5F7O1SteadyState.h"
// #include "DwarfElephantRBStructuresT6F1O1SteadyState.h"
// #include "DwarfElephantRBStructuresT7F8O80SteadyState.h"
// #include "DwarfElephantRBStructuresT7F8O2400SteadyState.h"
// #include "DwarfElephantRBStructuresT7F8O4747SteadyState.h"
// #include "DwarfElephantRBStructuresT8F9O1SteadyState.h"
// #include "DwarfElephantRBStructuresT8F9O10SteadyState.h"
// #include "DwarfElephantRBStructuresT8F9O20SteadyState.h"
// #include "DwarfElephantRBStructuresT8F9O980SteadyState.h"
// #include "DwarfElephantRBStructuresT8F9O2347SteadyState.h"
// #include "DwarfElephantRBStructuresT9F2O80SteadyState.h"
// #include "DwarfElephantRBStructuresT9F10O80SteadyState.h"
// #include "DwarfElephantRBStructuresT9F10O900SteadyState.h"
// #include "DwarfElephantRBStructuresT9F10O980SteadyState.h"
// #include "DwarfElephantRBStructuresT11F25O80SteadyState.h"
// #include "DwarfElephantRBStructuresT11F25O980SteadyState.h"
// #include "DwarfElephantRBStructuresT12F13O1SteadyState.h"
// #include "DwarfElephantRBStructuresT12F13O20SteadyState.h"
// #include "DwarfElephantRBStructuresT12F13O80SteadyState.h"
// #include "DwarfElephantRBStructuresT12F13O84SteadyState.h"
// #include "DwarfElephantRBStructuresT12F13O984SteadyState.h"
// #include "DwarfElephantRBStructuresT12F13O2347SteadyState.h"
// #include "DwarfElephantRBStructuresT12F28O980SteadyState.h"
// #include "DwarfElephantRBStructuresT13F13O84SteadyState.h"
// #include "DwarfElephantRBStructuresT14F14O30SteadyState.h"
// #include "DwarfElephantRBStructuresT14F14O84SteadyState.h"
// #include "DwarfElephantRBStructuresT14F14O983SteadyState.h"
// #include "DwarfElephantRBStructuresT15F8O80SteadyState.h"
// #include "DwarfElephantRBStructuresT15F16O80SteadyState.h"
// #include "DwarfElephantRBStructuresT16F16O80SteadyState.h"
#include "DwarfElephantRBT2F2O0SteadyStateExpansion.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;

  class EquationSystems;
  class RBConstruction;
  class RBEvaluation;
}

class DwarfElephantInitializeRBSystemSteadyState;

///In this class the subclasse of RBConstruction class is introduced.
class DwarfElephantRBConstructionSteadyState : public RBConstruction
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBConstructionSteadyState (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in);

  // Destructor
  virtual ~DwarfElephantRBConstructionSteadyState () { }

  // Type of the system
  typedef DwarfElephantRBConstructionSteadyState _sys_type;

  // Type of the parent
  typedef RBConstruction Parent;

  // Initialize data structure
  virtual void init_data();

  Real custom_train_reduced_basis(const bool resize_rb_eval_data = true);

  Real compute_residual_dual_norm(const unsigned int N);

  Real truth_solve(int plot_solution);

  unsigned int u_var;
  unsigned int lm_var;

};

///In this class the subclasse of RBEvaluation class is introduced. NOTE: ENSURE THAT THE CLASS IS USING THE CORRECT RBSTRUCTURES.
class DwarfElephantRBEvaluationSteadyState : public RBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluationSteadyState(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem);

    virtual Real get_stability_lower_bound();

// Outcommented at the moment, please remove comment marks in case you want to use the slower but less error
// prone error bound. Please, keep in mind that the slower error bound works in serial only.

//  Real rb_solve(unsigned int N)
//{
//  LOG_SCOPE("rb_solve()", "RBEvaluation");
//
//  if(N > get_n_basis_functions())
//    libmesh_error_msg("ERROR: N cannot be larger than the number of basis functions in rb_solve");
//
//  const RBParameters & mu = get_parameters();
//
//  // Resize (and clear) the solution vector
//  RB_solution.resize(N);
//
//  // Assemble the RB system
//  DenseMatrix<Number> RB_system_matrix(N,N);
//  RB_system_matrix.zero();
//
//  DenseMatrix<Number> RB_Aq_a;
//  for(unsigned int q_a=0; q_a<_rb_theta_expansion.get_n_A_terms(); q_a++)
//    {
//      RB_Aq_vector[q_a].get_principal_submatrix(N, RB_Aq_a);
//
//      RB_system_matrix.add(_rb_theta_expansion.eval_A_theta(q_a, mu), RB_Aq_a);
//    }
//
//  // Assemble the RB rhs
//  DenseVector<Number> RB_rhs(N);
//  RB_rhs.zero();
//
//  DenseVector<Number> RB_Fq_f;
//  for(unsigned int q_f=0; q_f<_rb_theta_expansion.get_n_F_terms(); q_f++)
//    {
//      RB_Fq_vector[q_f].get_principal_subvector(N, RB_Fq_f);
//      RB_rhs.add(_rb_theta_expansion.eval_F_theta(q_f, mu), RB_Fq_f);
//    }
//
//  // Solve the linear system
//  if(N > 0)
//    {
//      RB_system_matrix.lu_solve(RB_rhs, RB_solution);
//    }
//
//  // Evaluate RB outputs
//  DenseVector<Number> RB_output_vector_N;
//  for(unsigned int n=0; n<_rb_theta_expansion.get_n_outputs(); n++)
//    {
//      RB_outputs[n] = 0.;
//      for(unsigned int q_l=0; q_l<_rb_theta_expansion.get_n_output_terms(n); q_l++)
//        {
//          RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
//          RB_outputs[n] += _rb_theta_expansion.eval_output_theta(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
//        }
//    }
//
//  if(evaluate_RB_error_bound) // Calculate the error bounds
//    {
//      DwarfElephantRBConstructionSteadyState & sys_rb = fe_problem.es().get_system<DwarfElephantRBConstructionSteadyState>("RBSystem");
//      // Evaluate the dual norm of the residual for RB_solution_vector
//
////      // slower but less error prone error bound (does not work in parallel)
////      sys_rb.compute_residual_dual_norm(N);
//
//      // faster but more error prone error bound (does work in parallel)
//      compute_residual_dual_norm(N);
//
//      // Get lower bound for coercivity constant
//      const Real alpha_LB = get_stability_lower_bound();
//      // alpha_LB needs to be positive to get a valid error bound
//      libmesh_assert_greater ( alpha_LB, 0. );
//
//      // Evaluate the (absolute) error bound
//      Real abs_error_bound = epsilon_N / residual_scaling_denom(alpha_LB);
//
//      // Now compute the output error bounds
//      for(unsigned int n=0; n<_rb_theta_expansion.get_n_outputs(); n++)
//        {
//          RB_output_error_bounds[n] = abs_error_bound * eval_output_dual_norm(n, mu);
//        }
//
//      return abs_error_bound;
//    }
//  else // Don't calculate the error bounds
//    {
//      // Just return -1. if we did not compute the error bound
//      return -1.;
//    }
//}

  FEProblemBase & get_fe_problem(){return fe_problem;}

  FEProblemBase & fe_problem;
  // DwarfElephantRBT4F1O32SteadyStateExpansion _rb_theta_expansion;
  DwarfElephantRBT2F2O0SteadyStateExpansion _rb_theta_expansion;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSESSTEADYSTATE_H
