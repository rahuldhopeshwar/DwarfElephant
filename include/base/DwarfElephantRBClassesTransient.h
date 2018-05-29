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
#include "FEProblemBase.h"
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemTransient.h"
#include "DwarfElephantOfflineOnlineStageTransient.h"

#include "DwarfElephantRBStructuresT1F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT2F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT3F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT3F1O3M1Transient.h"
#include "DwarfElephantRBStructuresT3F4O1M2Transient.h"
#include "DwarfElephantRBStructuresT3F4O3M2Transient.h"
#include "DwarfElephantRBStructuresT4F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT5F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT5F4O1M2Transient.h"

// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;

  class EquationSystems;
  class TransientRBConstruction;
  class TransientRBEvaluation;
}

///-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
class DwarfElephantRBConstructionTransient : public TransientRBConstruction
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBConstructionTransient (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in);

  // Destructor
  virtual ~DwarfElephantRBConstructionTransient () { }

  // Type of the system
  typedef DwarfElephantRBConstructionTransient _sys_type;

  // Type of the parent
  typedef TransientRBConstruction Parent;

  // Initialize data structure
  virtual void init_data() override;

  // virtual Real truth_solve(int write_interval) libmesh_override;

  virtual void print_info() override;
  //
  // virtual void initialize_truth() override;

  // void add_IC_to_RB_space();

  void update_RB_initial_condition_all_N();

  // virtual void update_system() override;

  virtual Real get_RB_error_bound() override;

  virtual Real train_reduced_basis(const bool resize_rb_eval_data=true) override;
  virtual Real train_reduced_basis_steady(const bool resize_rb_eval_data=true);

  unsigned int u_var;

};

///------------------------DWARFELEPHANTRBEVALUATION------------------------
class DwarfElephantRBEvaluationTransient : public TransientRBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluationTransient(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem);

  virtual Real get_stability_lower_bound() override;

  FEProblemBase & get_fe_problem() {return fe_problem;}

  FEProblemBase & fe_problem;
  DwarfElephantRBT3F4O1M2TransientExpansion _rb_theta_expansion;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSESTRANSIENT_H
