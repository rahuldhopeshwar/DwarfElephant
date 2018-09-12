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

#include <fstream>
// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"

// libMesh includes (RB package)
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_eim_assembly.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// MOOSE includes
#include "FEProblemBase.h"

// MOOSE includes (DwarfElephant package)

#include "DwarfElephantRBStructuresT1F1O1SteadyState.h"
#include "DwarfElephantRBStructuresT2F1O1SteadyState.h"
#include "DwarfElephantRBStructuresT3F1O1SteadyState.h"
#include "DwarfElephantRBStructuresT3F3O1SteadyState.h"
#include "DwarfElephantRBStructuresT3F1O3SteadyState.h"
#include "DwarfElephantRBStructuresT4F1O1SteadyState.h"
#include "DwarfElephantRBStructuresT5F1O1SteadyState.h"
#include "DwarfElephantRBStructuresT5F3O1SteadyState.h"
#include "DwarfElephantRBStructuresT6F1O1SteadyState.h"
#include "DwarfElephantEIMStructures.h"
///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class NumericVector;
  template <typename T> class PetscMatrix;

  class EquationSystems;
  class RBConstruction;
  class RBEvaluation;
  class RBEIMConstruction;
  class RBEIMEvaluation;
  class RBEIMAssembly;
}

class DwarfElephantInitializeRBSystemSteadyState;

///-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
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
  NumericVector<Number> * get_nonAffineF() // To test against EIM example from Martin's publication
  {
    return _nonAffineF.get();
  }

  SparseMatrix<Number> * get_nonAffineA() // To test against EIM example from Martin's publication
  {
    return _nonAffineA.get();
  }
//  Real train_reduced_basis(const bool resize_rb_eval_data = true);

  void allocate_EIM_error_structures() // To test against EIM example from Martin's publication
  {
    _nonAffineA = SparseMatrix<Number>::build(this->comm());
    libMesh::DofMap & dof_map = this -> get_dof_map();
    dof_map.attach_matrix(*_nonAffineA);
    _nonAffineA->init();
    _nonAffineA->zero();

    _nonAffineF = NumericVector<Number>::build(this->comm());
    _nonAffineF->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  }

  Real compute_residual_dual_norm(const unsigned int N);
  
  void TestEIMAccuracy()// To test against EIM example from Martin's publication
  {
    RBParameters mu;
    mu.set_value("mu_0", -0.01);
    mu.set_value("mu_1", -0.01);
    std::ofstream EIMErrorFile;
    EIMErrorFile.open("EIMError_vs_M.csv");
   
    EIMErrorFile << "M, EIMError" << std::endl;
    for (unsigned int M = 1; M < get_rb_theta_expansion().get_n_F_terms(); M++)
      EIMErrorFile << M << ", " << find_EIMFE_error(M,mu) << std::endl;

    EIMErrorFile.close();
  }
  Real find_EIMFE_error(unsigned int M, RBParameters mu)// To test against EIM example from Martin's publication
  {
    // compute full FE solution
      // assemble full FE system
      // A_matrix = Aq[0] + _nonAffineA
      // F_vector = _nonAffineF
      std::unique_ptr<NumericVector<Number>> fullFEsolution, EIM_FEsolution;
      fullFEsolution -> zero();
      fullFEsolution -> close();

      EIM_FEsolution -> zero();
      EIM_FEsolution -> close();


      this -> matrix -> zero();
      this -> rhs -> zero();
      
      this -> matrix -> close();
      this -> rhs -> close();
      
      //matrix -> add(1.0,*get_Aq(0));
      matrix -> add(1.0,*_nonAffineA);
      rhs -> add(*_nonAffineF);

      this -> matrix -> close();
      this -> rhs -> close();
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
    fullFEsolution -> add(*solution);
      



     // solve equation

    // compute EIM_FE solution for 1<=M<=M_max
      // assemble EIM_FE system for M
      // solve equation


      this -> matrix -> zero();
      this -> rhs -> zero();
      
      this -> matrix -> close();
      this -> rhs -> close();
      
      for (unsigned int q_a=0; q_a<M+1; q_a++)
      {
        matrix->add(get_rb_theta_expansion().eval_A_theta(q_a, mu), *get_Aq(q_a));
      }

    std::unique_ptr<NumericVector<Number>> temp_vec = NumericVector<Number>::build(this->comm());
    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    for (unsigned int q_f=0; q_f<M; q_f++)
      {
        *temp_vec = *get_Fq(q_f);
        temp_vec->scale( get_rb_theta_expansion().eval_F_theta(q_f, mu) );
        rhs->add(*temp_vec);
      }

      this -> matrix -> close();
      this -> rhs -> close();

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
    EIM_FEsolution -> add(*solution);
    *fullFEsolution -= *EIM_FEsolution;
    get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector, *fullFEsolution);
    Number truthError_X_norm = std::sqrt(inner_product_storage_vector->dot(*fullFEsolution));

    return truthError_X_norm;
   // Compute X inner product of the difference between full FE solution and EIM_FE solution
  }
  

  unsigned int u_var;
  std::unique_ptr<SparseMatrix<Number>> _nonAffineA; // To test against EIM example from Martin's publication
  std::unique_ptr<NumericVector<Number>> _nonAffineF; // To test against EIM example from Martin's publication

};

///------------------------DWARFELEPHANTRBEVALUATION------------------------
class DwarfElephantRBEvaluationSteadyState : public RBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluationSteadyState(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem);

  virtual ~DwarfElephantRBEvaluationSteadyState() {}
    virtual Real get_stability_lower_bound();

  FEProblemBase & get_fe_problem(){return fe_problem;}

  FEProblemBase & fe_problem;
  DwarfElephantEIMTestRBThetaExpansion _eim_test_rb_theta_expansion;
};


class DwarfElephantEIMEvaluationSteadyState : public RBEIMEvaluation
{
public:

  DwarfElephantEIMEvaluationSteadyState(const libMesh::Parallel::Communicator & comm);
  
  virtual ~DwarfElephantEIMEvaluationSteadyState() {}

  ShiftedGaussian sg;
};

// A simple subclass of RBEIMConstruction.
class DwarfElephantEIMConstructionSteadyState : public RBEIMConstruction
{
public:

  /**
   * Constructor.
   */
  DwarfElephantEIMConstructionSteadyState (EquationSystems & es,
                         const std::string & name_in,
                         const unsigned int number_in);

  
  virtual ~DwarfElephantEIMConstructionSteadyState() {}
  /**
   * The type of the parent.
   */
  typedef RBEIMConstruction Parent;

  virtual std::unique_ptr<ElemAssembly> build_eim_assembly(unsigned int index);
  
  virtual void init_data();
  
  /**
   * Initialize the implicit system that is used to perform L2 projections.
   */
  virtual void init_implicit_system();

  /**
   * Initialize the explicit system that is used to store the basis functions.
   */
  virtual void init_explicit_system();

  virtual Real train_reduced_basis(const bool resize_rb_eval_data=true)
{
  // precompute all the parametrized functions that we'll use in the greedy
  initialize_parametrized_functions_in_training_set();

  LOG_SCOPE("train_reduced_basis()", "RBConstruction");

  int count = 0;

  // initialize rb_eval's parameters
  get_rb_evaluation().initialize_parameters(*this);

  // possibly resize data structures according to Nmax
  if (resize_rb_eval_data)
    {
      get_rb_evaluation().resize_data_structures(get_Nmax());
    }

  // Clear the Greedy param list
  for (std::size_t i=0; i<get_rb_evaluation().greedy_param_list.size(); i++)
    get_rb_evaluation().greedy_param_list[i].clear();

  get_rb_evaluation().greedy_param_list.clear();

  Real training_greedy_error;


  // If we are continuing from a previous training run,
  // we might already be at the max number of basis functions.
  // If so, we can just return.
  if (get_rb_evaluation().get_n_basis_functions() >= get_Nmax())
    {
      libMesh::out << "Maximum number of basis functions reached: Nmax = "
                   << get_Nmax() << std::endl;
      return 0.;
    }


  // Compute the dual norms of the outputs if we haven't already done so
  compute_output_dual_innerprods();
std::ofstream greedyFile;
greedyFile.open("EIM_Example_Martin_InterpolationPoints_L2bestfit.csv");
  // Compute the Fq Riesz representor dual norms if we haven't already done so
  compute_Fq_representor_innerprods();

  libMesh::out << std::endl << "---- Performing Greedy basis enrichment ----" << std::endl;
  Real initial_greedy_error = 0.;
  bool initial_greedy_error_initialized = false;
  while (true)
    {
      libMesh::out << std::endl << "---- Basis dimension: "
                   << get_rb_evaluation().get_n_basis_functions() << " ----" << std::endl;

      if (count > 0 || (count==0 && use_empty_rb_solve_in_greedy))
        {
          libMesh::out << "Performing RB solves on training set" << std::endl;
          training_greedy_error = compute_max_error_bound();
          
          
          libMesh::out << "Maximum error bound is " << training_greedy_error << std::endl << std::endl;

          // record the initial error
          if (!initial_greedy_error_initialized)
            {
              initial_greedy_error = training_greedy_error;
              initial_greedy_error_initialized = true;
            }

          // Break out of training phase if we have reached Nmax
          // or if the training_tolerance is satisfied.
          if (greedy_termination_test(training_greedy_error, initial_greedy_error, count))
            break;
        }

      libMesh::out << "Performing truth solve at parameter:" << std::endl;
      print_parameters();
      
      
      // Update the list of Greedily selected parameters
      this->update_greedy_param_list();

      // Perform an Offline truth solve for the current parameter
      truth_solve(-1);
      
      // Add orthogonal part of the snapshot to the RB space
      libMesh::out << "Enriching the RB space" << std::endl;
      enrich_RB_space();
      
      update_system();

      // Increment counter
      count++;
    }
   RBEIMEvaluation & eim_eval = cast_ref<RBEIMEvaluation &>(get_rb_evaluation());
   for (unsigned int i = 0; i < count ; i++) { eim_eval.interpolation_points[i].write_unformatted(greedyFile); }
   greedyFile.close();
  this->update_greedy_param_list();

  return training_greedy_error;

}
  /**
   * Variable number for u.
   */
  unsigned int u_var;

  std::vector<std::unique_ptr<DwarfElephantEIMFAssembly>> _rb_eim_assembly_objects_new;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSESSTEADYSTATE_H
