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
#include "DwarfElephantRBClasses.h"


///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;
  template <typename T> class PetscVector;
}

class MooseMesh;
class DwarfElephantRBConstruction;
class DwarfElephantInitializeRBSystem;

template<>
InputParameters validParams<DwarfElephantInitializeRBSystem>();

class DwarfElephantInitializeRBSystem :
  public GeneralUserObject
{
  public:
    DwarfElephantInitializeRBSystem(const InputParameters & params);

    void onlineStage();

    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;

    DwarfElephantRBConstruction * _rb_con_ptr;
    std::vector<SparseMatrix<Number> *> _jacobian_subdomain;
    std::vector<NumericVector<Number> *> _residuals;
    std::vector<NumericVector<Number> *> _outputs;

  protected:
    bool _use_displaced;
    bool _skip_matrix_assembly_in_rb_system;
    bool _skip_vector_assembly_in_rb_system;
    bool _offline_stage;
    bool _online_stage;
    bool _F_equal_to_output;
    bool _store_basis_functions;

    unsigned int _online_N;
    unsigned int _qa;
    unsigned int _qf;
    unsigned int _ql;

    Real _online_mu;

    std::string _parameters_filename;
    std::string _system_name;

    EquationSystems & _es;
    TransientNonlinearImplicitSystem & _sys;

    MooseMesh * _mesh_ptr;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTINITIALIZERBSYSTEM_H
