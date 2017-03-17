///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTPREPARERBSYSTEM_H
#define DWARFELEPHANTPREPARERBSYSTEM_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"

// MOOSE includes
#include "NodalUserObject.h"
#include "DisplacedProblem.h"
#include "Assembly.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;
}

class Assembly;
class DisplacedProblem;
class DwarfElephantPrepareRBSystem;

template<>
InputParameters validParams<DwarfElephantPrepareRBSystem>();

class DwarfElephantPrepareRBSystem :
  public NodalUserObject
{
  public:
    DwarfElephantPrepareRBSystem(const InputParameters & params);

    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;
    virtual void threadJoin(const UserObject & y);


  protected:
    DofMap & _dof_map;
    bool _use_displaced;
    EquationSystems & _es;

//    UniquePtr<SparseMatrix <Number> > _Aq_qa;
    SparseMatrix <Number> * _Aq_qa;

    TransientNonlinearImplicitSystem & _sys;

};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTPREPARERBSYSTEM_H
