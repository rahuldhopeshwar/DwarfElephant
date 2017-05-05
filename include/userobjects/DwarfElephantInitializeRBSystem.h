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
#include "DwarfElephantRBClassesSteadyState.h"
#include "DwarfElephantRBClassesTransient.h"
#include "CacheBoundaries.h"


///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  class RBConstruction;
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;
  template <typename T> class PetscVector;
}

class MooseMesh;
class DwarfElephantRBConstruction;
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
    void initVariable();

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

    std::string _parameters_filename;
    std::string _rb_variable_name;

    EquationSystems & _es;
    MooseMesh * _mesh_ptr;
    DwarfElephantRBConstructionSteadyState * _rb_con_ptr;

    SparseMatrix <Number> * _inner_product_matrix;
    std::vector<SparseMatrix <Number> *> _jacobian_subdomain;
    std::vector<NumericVector <Number> *> _residuals;
    std::vector<NumericVector <Number> *> _outputs;

    std::vector <numeric_index_type> _cached_jacobian_subdomain_contribution_rows;
    std::vector <numeric_index_type> _cached_jacobian_subdomain_contribution_cols;
    std::vector <Real> _cached_jacobian_subdomain_contribution_vals;

    const std::vector<ExecFlagType> & _exec_flags;

    friend class RBKernel;
    friend class RBNodalBC;
    friend class DwarfElephantOfflineOnlineStage;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTINITIALIZERBSYSTEM_H
