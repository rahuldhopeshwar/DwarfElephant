///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTINITIALIZERBSYSTEM_H
#define DWARFELEPHANTINITIALIZERBSYSTEM_H

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
#include "CacheBoundaries.h"


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
class DwarfElephantInitializeRBSystem;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantInitializeRBSystem>();

///-------------------------------------------------------------------------
class DwarfElephantInitializeRBSystem :
  public GeneralUserObject
{

//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantInitializeRBSystem(const InputParameters & params);

    /* Methods */
    void initializeOfflineStage();

    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;

//--------------------------------PROTECTED---------------------------------
  protected:

    /* Attributes */
    bool _use_displaced;
    bool _skip_matrix_assembly_in_rb_system;
    bool _skip_vector_assembly_in_rb_system;
    bool _offline_stage;
    bool _compliant;

    unsigned int _qa;
    unsigned int _qf;
    unsigned int _ql;

    std::string _system_name;
    std::string _parameters_filename;

    EquationSystems & _es;
    MooseMesh * _mesh_ptr;
    TransientNonlinearImplicitSystem * _sys;
    DwarfElephantRBConstructionSteadyState * _rb_con_ptr;

    SparseMatrix <Number> * _inner_product_matrix;
    std::vector<SparseMatrix <Number> *> _jacobian_subdomain;
    std::vector<NumericVector <Number> *> _residuals;
    std::vector<NumericVector <Number> *> _outputs;

    const std::vector<ExecFlagType> & _exec_flags;

    Function * _function;
    CacheBoundaries * _cache_boundaries;


    friend class RBKernel;
    friend class RBNodalBC;
    friend class DwarfElephantOfflineOnlineStage;
    friend class DwarfElephantRBConstructionSteadyState;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTINITIALIZERBSYSTEM_H
