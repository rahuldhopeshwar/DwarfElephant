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
#include "MooseMesh.h"
#include "Assembly.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
  template <typename T> class PetscMatrix;
}

class MooseMesh;
class Assembly;
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
    bool _use_displaced;

    std::string _system_name;

    const std::set<SubdomainID> & _block_ids;

    EquationSystems & _es;
    TransientNonlinearImplicitSystem & _sys;
    DofMap & _dof_map;

    MooseMesh * _mesh_ptr;

    std::map<SubdomainID, UniquePtr< NumericVector<Number> > > _Fq_a;
//    std::map<SubdomainID, UniquePtr< SparseMatrix<Number> > > _Aq_a;
    std::map<SubdomainID, SparseMatrix<Number> * > _Aq_a;
    UniquePtr<SparseMatrix<Number>> _A_test;

};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTPREPARERBSYSTEM_H
