///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTOFFLINESTAGE_H
#define DWARFELEPHANTOFFLINESTAGE_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/getpot.h"

// MOOSE includes
#include "NodalUserObject.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBClasses.h"
#include "DwarfElephantRBClassesAssemble.h"


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
class DwarfElephantOfflineStage;

template<>
InputParameters validParams<DwarfElephantOfflineStage>();

class DwarfElephantOfflineStage :
  public NodalUserObject // Code does not work for GeneralUserObject
{
  public:
    DwarfElephantOfflineStage(const InputParameters & params);

    void offlineStage();
    void transferAffineOperators(bool _skip_matrix_assembly_in_rb_system, bool _skip_vector_assembly_in_rb_system);

    virtual void initialize() override;
    virtual void execute() override;
    virtual void threadJoin(const UserObject & y);
    virtual void finalize() override;

  protected:
    bool _use_displaced;
    bool _skip_matrix_assembly_in_rb_system;
    bool _skip_vector_assembly_in_rb_system;
    bool _F_equal_to_output;
    bool _store_basis_functions;

    std::string _parameters_filename;
    std::string _system_name;

    const std::set<SubdomainID> & _block_ids;

    unsigned int _qa;
    unsigned int _qf;
    unsigned int _ql;

    EquationSystems & _es;
    TransientNonlinearImplicitSystem & _sys;

    MooseMesh * _mesh_ptr;
    DwarfElephantRBConstruction * _rb_con_ptr;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTINITIALIZERBSYSTEM_H
