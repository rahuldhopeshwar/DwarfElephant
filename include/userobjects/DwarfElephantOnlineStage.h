///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTONLINESTAGE_H
#define DWARFELEPHANTONLINESTAGE_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"

// MOOSE includes
#include "NodalUserObject.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"
#include "NodalBC.h"

#include "DwarfElephantInitializeRBSystem.h"

///-------------------------------------------------------------------------
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
}
// Forward Declarations
class DisplacedProblem;
class MooseMesh;
class DwarfElephantOnlineStage;

template<>
InputParameters validParams<DwarfElephantOnlineStage>();

class DwarfElephantOnlineStage :
  public NodalUserObject
{
  public:
    DwarfElephantOnlineStage(const InputParameters & params);

    void onlineStage();

    virtual void initialize() override;
    virtual void execute() override;
    virtual void threadJoin(const UserObject & y) override;
    virtual void finalize() override;

  protected:
    bool _use_displaced;
    std::string _system_name;
    EquationSystems & _es;
    TransientNonlinearImplicitSystem & _sys;
    MooseMesh * _mesh_ptr;
    MooseObjectWarehouse <NodalBC> _nodal_bcs;

    NonlinearVariableName _variable_name;
    MooseVariable & _var;
   const DwarfElephantInitializeRBSystem & _initialize_rb_system;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTONLINESTAGE_H
