///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTINITIALIZERBSYSTEMSTEADYSTATE_H
#define DWARFELEPHANTINITIALIZERBSYSTEMSTEADYSTATE_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/getpot.h"

// MOOSE includes
#include "GeneralUserObject.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantSystem.h"
#include "DwarfElephantRBClassesSteadyState.h"
//#include "DwarfElephantRBConstructionSteadyState.h"
//#include "DwarfElephantRBClassesTransient.h"


///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
//  class RBConstruction;
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;
  template <typename T> class PetscVector;
}

class MooseMesh;
class DwarfElephantSystem;
class DwarfElephantRBConstructionSteadyState;
class DwarfElephantInitializeRBSystemSteadyState;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantInitializeRBSystemSteadyState>();

///-------------------------------------------------------------------------
class DwarfElephantInitializeRBSystemSteadyState :
  public GeneralUserObject
{

//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantInitializeRBSystemSteadyState(const InputParameters & params);

    /* Methods */

    // Initializes all required matrices and vectors for the RB solve.
    void initializeOfflineStage();

    // Initializes the RB System.
    virtual void initialize() override;

    // Method not used in this UserObject.
    virtual void execute() override;

    // Method not used in this UserObject.
    virtual void finalize() override;

//--------------------------------PROTECTED---------------------------------
  protected:

    /* Attributes */
    bool _use_displaced;
    bool _skip_matrix_assembly_in_rb_system;
    bool _skip_vector_assembly_in_rb_system;
    bool _offline_stage;
    bool _compliant;

    unsigned int _n_outputs;
    unsigned int _qa;
    unsigned int _qf;
    std::vector<unsigned int> _ql;

    std::string _system_name;
    std::string _parameters_filename;

    EquationSystems & _es;
    MooseMesh * _mesh_ptr;
    TransientNonlinearImplicitSystem * _sys;
    DwarfElephantRBConstructionSteadyState * _rb_con_ptr;

    SparseMatrix <Number> * _inner_product_matrix;
    std::vector<SparseMatrix <Number> *> _jacobian_subdomain;
    std::vector<NumericVector <Number> *> _residuals;
    std::vector<std::vector<NumericVector <Number> *> > _outputs;

    const std::vector<ExecFlagType> & _exec_flags;


    friend class RBKernel;
    friend class DwarfElephantRBNodalBC;
    friend class DwarfElephantRBIntegratedBC;
    friend class DwarfElephantOfflineOnlineStageSteadyState;
    friend class DwarfElephantRBEvaluationSteadyState;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTINITIALIZERBSYSTEMSTEADYSTATE_H
