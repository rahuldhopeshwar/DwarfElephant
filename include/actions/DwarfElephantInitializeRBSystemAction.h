///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTINITIALIZERBSYSTEMACTION_H
#define DWARFELEPHANTINITIALIERBSYSTEMACTION_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/getpot.h"

// MOOSE includes
#include "Action.h"
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
class DwarfElephantInitializeRBSystemAction;

template<>
InputParameters validParams<DwarfElephantInitializeRBSystemAction>();

class DwarfElephantInitializeRBSystemAction :
  public Action
{
  public:
    DwarfElephantInitializeRBSystemAction(InputParameters params);

    void initializeRBSystem();

    virtual void act() override;

  protected:
    bool _use_displaced;
    bool _skip_matrix_assembly_in_rb_system;
    bool _skip_vector_assembly_in_rb_system;

    std::string _parameters_filename;

    EquationSystems & _es;

    MooseMesh * _mesh_ptr;
    DwarfElephantRBConstruction * _rb_con_ptr;

};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTINITIALIZERBSYSTEMACTION_H
