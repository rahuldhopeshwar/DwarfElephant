/**
 * In this class simplified subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
 *
 * RBSimpleEvaluation: requires only the definition of the lower coercivity
 * constant. The value is here specified for a Conduction problem.
 *
 * RBSimpleConstruction: In order to construct the RB System with the
 * RBSimpleEvaluation subclass the method build_rb_evaluation needs to be
 * overriden.
 *
 * NOTE: ENSURE THAT THE CLASS IS INHERITING FROM THE CORRECT RBSTRUCTURES.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBCLASSES_H
#define DWARFELEPHANTRBCLASSES_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/fe_base.h"

//libMesh includes (RB package)
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"

#include "FEProblemBase.h"

//MOOSE includes (DwarfElephant package)
#include "RBStructuresP1_3ThetaEqualMuCompliant.h"

using libMesh::DirichletBoundary;
using libMesh::EquationSystems;
using libMesh::FEMContext;
using libMesh::RBConstruction;
using libMesh::RBEvaluation;
using libMesh::Real;
using libMesh::UniquePtr;

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

  virtual Real get_stability_lower_bound(){ return 0.05; };

  RBP1_3ThetaEqualMuThetaCompliantExpansion _rb_theta_expansion;
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
 virtual void init_data()
  {
    u_var = this->add_variable ("u", libMesh::FIRST);

    // Generate a DirichletBoundary object
    dirichlet_bc = build_zero_dirichlet_boundary_object();

    // Set the Dirichet boundary IDs
    // and the Dirichlet boundary variable numbers
//    dirichlet_bc->b.insert(0);
    dirichlet_bc->b.insert(1);
//    dirichlet_bc->b.insert(2);
    dirichlet_bc->b.insert(3);
    dirichlet_bc->variables.push_back(u_var);

    // Attach dirichlet_bc (must do this _before_ Parent::init_data)
    get_dof_map().add_dirichlet_boundary(*dirichlet_bc);

    Parent::init_data();

    // Set the rb_assembly_expansion for this Construction object.
    // The theta expansion comes from the RBEvaluation object.
    set_rb_assembly_expansion(_rb_assembly_expansion);

    // We need to define an inner product matrix for this problem
    set_inner_product_assembly(_rb_assembly_expansion.A0_assembly);
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

  /**
   * Variable number for u.
   */
  unsigned int u_var;

  /**
   * The object that stores the "assembly" expansion of the parameter dependent PDE,
   * i.e. the objects that define how to assemble the set of parameter-independent
   * operators in the affine expansion of the PDE.
   */
  RBP1_3ThetaEqualMuAssemblyCompliantExpansion _rb_assembly_expansion;

  /**
   * The object that defines which degrees of freedom are on a Dirichlet boundary.
   */
  UniquePtr<DirichletBoundary> dirichlet_bc;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSES_H

//Real compute_max_error_bound()
//{
//  LOG_SCOPE("compute_max_error_bound()", "DwarfElephantRBConstruction");
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
//Real train_reduced_basis(const bool resize_rb_eval_data=true)
//{
//  LOG_SCOPE("train_reduced_basis()", "DwarfElephantRBConstruction");
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
//  for(unsigned int i=0; i<get_rb_evaluation().greedy_param_list.size(); i++)
//    {
//      get_rb_evaluation().greedy_param_list[i].clear();
//    }
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
/////__________________________________________________________________________________________
///// Use Moose Postprocessor VariableInnerProduct ?
////  // Compute the dual norms of the outputs if we haven't already done so
////  compute_output_dual_innerprods();
//
////  // Compute the Fq Riesz representor dual norms if we haven't already done so
////  compute_Fq_representor_innerprods();
/////__________________________________________________________________________________________
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
////
////          // record the initial error
////          if (!initial_greedy_error_initialized)
////            {
////              initial_greedy_error = training_greedy_error;
////              initial_greedy_error_initialized = true;
////            }
////
////          // Break out of training phase if we have reached Nmax
////          // or if the training_tolerance is satisfied.
////          if (greedy_termination_test(training_greedy_error, initial_greedy_error, count))
//            break;
//        }
////
////      libMesh::out << "Performing truth solve at parameter:" << std::endl;
////      print_parameters();
////
////      // Update the list of Greedily selected parameters
////      this->update_greedy_param_list();
////
////      // Perform an Offline truth solve for the current parameter
////      truth_solve(-1);
////
////      // Add orthogonal part of the snapshot to the RB space
////      libMesh::out << "Enriching the RB space" << std::endl;
////      enrich_RB_space();
////
////      update_system();
////
//      // Increment counter
//      count++;
//    }
////  this->update_greedy_param_list();
////
//  return training_greedy_error;
//}
