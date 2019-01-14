/**
 * In this class subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
 *
 * NOTE: ENSURE THAT THE CLASS IS USING THE CORRECT RBSTRUCTURES.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBCLASSES4DVAR_H
#define DWARFELEPHANTRBCLASSES4DVAR_H

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
#include "DwarfElephantRBClassesTransient.h"

#include "DwarfElephantRBStructuresT1F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT2F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT2F2O1M2Transient.h"
#include "DwarfElephantRBStructuresT2F3O1M2Transient.h"
#include "DwarfElephantRBStructuresT2F3O3M2Transient.h"
#include "DwarfElephantRBStructuresT3F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT3F1O3M1Transient.h"
#include "DwarfElephantRBStructuresT3F1O80M1Transient.h"
#include "DwarfElephantRBStructuresT3F4O1M2Transient.h"
#include "DwarfElephantRBStructuresT3F4O3M2Transient.h"
#include "DwarfElephantRBStructuresT3F3O3M2IC1Transient.h"
#include "DwarfElephantRBStructuresT4F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT4F1O1M1IC1Transient.h"
#include "DwarfElephantRBStructuresT5F1O1M1Transient.h"
#include "DwarfElephantRBStructuresT5F4O1M2Transient.h"
#include "DwarfElephantRBStructuresT6F1O1M1IC3Transient.h"

// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;

  class EquationSystems;
}

class DwarfElephantRBConstructionTransient;
class DwarfElephantRBEvaluationTransient;

///In this class the subclasse of TransientRBConstruction class is introduced.
class DwarfElephantRBConstruction4DVar : public DwarfElephantRBConstructionTransient
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBConstruction4DVar (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in);

  // Destructor
  virtual ~DwarfElephantRBConstruction4DVar () { }

  // Type of the system
  typedef DwarfElephantRBConstruction4DVar _sys_type;

  // Type of the parent
  typedef DwarfElephantRBConstructionTransient Parent;

  virtual void init_data() override;

  virtual Real train_reduced_basis(const bool resize_rb_eval_data=true) override;

  virtual void allocate_data_structures() override;

  // unsigned int get_n_qois() const {return n_qois;}
  //
  // virtual void set_n_qois(unsigned int n_qois_in);

  // std::vector<Real> get_qoi_weights() const {return qoi_weights;}

  // virtual void set_qoi_weights(std::vector<Real> qoi_weights_in);

protected:
  // unsigned int n_qois;
  // std::vector<Real> qoi_weights;
  std::vector<DenseVector<Number>> obs_data_all_k;

  // std::vector<std::unique_ptr<NumericVector<Number>>> Cq_vector;

friend class DwarfElephantInitializeRBSystem4DVar;

};

///In this class the subclasse of TransientRBEvaluation class is introduced. NOTE: ENSURE THAT THE CLASS IS USING THE CORRECT RBSTRUCTURES.
class DwarfElephantRBEvaluation4DVar : public DwarfElephantRBEvaluationTransient
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluation4DVar(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem);

  typedef DwarfElephantRBEvaluationTransient Parent;

  FEProblemBase & fe_problem;

  DwarfElephantRBT3F1O80M1TransientExpansion _rb_theta_expansion;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSES4DVAR_H
