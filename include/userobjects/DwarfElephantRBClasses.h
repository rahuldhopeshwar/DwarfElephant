/**
 * In this class simplified subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
 *
 * DwarfElephantRBEvaluation: requires only the definition of the lower
 * coercivity constant. The value is here specified for a Conduction
 * problem.
 *
 * DwarfElephantRBConstruction: In order to construct the RB System with the
 * RBSimpleEvaluation subclass the method build_rb_evaluation needs to be
 * overriden.
 *
 * NOTE: ENSURE THAT THE CLASS IS INHERITING FROM THE CORRECT RBSTRUCTURES.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBCLASSES_H
#define DWARFELEPHANTRBCLASSES_H

///---------------------------------INCLUDES--------------------------------
#include "libmesh/zero_function.h"

//libMesh includes (RB package)
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"

//MOOSE includes (DwarfElephant package)
#include "RBStructuresP1Theta3ThetaEqualMu.h"

// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  class RBConstruction;
  class RBEvaluation;
  class DirichletBoundary;
}

///---------------------------RBSIMPLEEVALUATION----------------------------
class DwarfElephantRBEvaluation : public RBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluation(const libMesh::Parallel::Communicator & comm):
    RBEvaluation(comm)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }

  virtual Real get_stability_lower_bound() { return 1;}

  RBP1Theta3ThetaEqualMuExpansion _rb_theta_expansion;
};

///--------------------------RBSIMPLECONSTRUCTION---------------------------
class DwarfElephantRBConstruction : public RBConstruction
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBConstruction (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : Parent(es, name_in, number_in),
      dirichlet_bc(UniquePtr<DirichletBoundary>())
  {}

  // Destructor
  virtual ~DwarfElephantRBConstruction () { }

  // Type of the system
  typedef DwarfElephantRBConstruction _sys_type;

  // Type of the parent
  typedef RBConstruction Parent;

  // Initialize data structure
  UniquePtr<DirichletBoundary> build_zero_dirichlet_boundary_object()
{
  ZeroFunction<> zf;
  ConstFunction<Number> one(1);

  std::set<boundary_id_type> dirichlet_ids;
  std::vector<unsigned int> variables;

  // The DirichletBoundary constructor clones zf, so it's OK that zf is only in local scope
  return UniquePtr<DirichletBoundary> (new DirichletBoundary(dirichlet_ids, variables, &zf));
}

 virtual void init_data()
  {
    u_var = this->add_variable ("u", libMesh::FIRST);

    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(2);

    std::vector<unsigned int> variables;
    variables.push_back(u_var);

    get_dof_map().add_dirichlet_boundary(DirichletBoundary(boundary_ids, variables, ConstFunction<Number>(1.)));

    Parent::init_data();

     set_rb_assembly_expansion(rb_assembly_expansion);

    // We need to define an inner product matrix for this problem
    set_inner_product_assembly(rb_assembly_expansion.A0_assembly);
  }

  /**
   * Pre-request all relevant element data.
   */
  virtual void init_context(FEMContext & c)
  {
    // For efficiency, we should prerequest all
    // the data we will need to build the
    // linear system before doing an element loop.
    FEBase * elem_fe = libmesh_nullptr;
    c.get_element_fe(u_var, elem_fe);

    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_dphi();
  }


//  Real train_reduced_basis(const bool resize_rb_eval_data)
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

  /**
   * Variable number for u.
   */
  unsigned int u_var;

  /**
   * The object that stores the "assembly" expansion of the parameter dependent PDE,
   * i.e. the objects that define how to assemble the set of parameter-independent
   * operators in the affine expansion of the PDE.
   */
  DwarfElephantRBAssemblyExpansion rb_assembly_expansion;

  /**
   * The object that defines which degrees of freedom are on a Dirichlet boundary.
   */

   UniquePtr<DirichletBoundary> dirichlet_bc;
};
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
//              if (!is_quiet())
////                libMesh::out << "Starting solve n=" << n << ", q_l=" << q_l
////                             << " in RBConstruction::compute_output_dual_innerprods() at "
////                             << Utility::get_timestamp() << std::endl;
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
//              if (!is_quiet())
//                {
////                  libMesh::out << "Finished solve n=" << n << ", q_l=" << q_l
////                               << " in RBConstruction::compute_output_dual_innerprods() at "
////                               << Utility::get_timestamp() << std::endl;
//
//                  libMesh::out << this->n_linear_iterations()
//                               << " iterations, final residual "
//                               << this->final_linear_residual() << std::endl;
//                }
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

//  unsigned int u_var;
//  UniquePtr<DirichletBoundary> dirichlet_bc;
//};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSES_H
