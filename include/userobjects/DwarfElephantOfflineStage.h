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
#include "GeneralUserObject.h"
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
  public GeneralUserObject
{
  public:
    DwarfElephantOfflineStage(const InputParameters & params);

    void setInnerProductMatrix();
    void offlineStage();
    void transferAffineOperators(bool _skip_matrix_assembly_in_rb_system, bool _skip_vector_assembly_in_rb_system);

    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;

  protected:
    bool _use_displaced;
    bool _store_basis_functions;
    bool _skip_matrix_assembly_in_rb_system;
    bool _skip_vector_assembly_in_rb_system;
    bool _compliant;

    std::string _parameters_filename;
    std::string _system_name;

    EquationSystems & _es;
    TransientNonlinearImplicitSystem & _sys;
    const DwarfElephantInitializeRBSystem & _initialize_rb_system;

    MooseMesh * _mesh_ptr;
    const std::set<SubdomainID> & _subdomain_ids;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTOFFLINESTAGE_H
