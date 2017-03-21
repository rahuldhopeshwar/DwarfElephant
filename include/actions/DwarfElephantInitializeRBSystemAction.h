#ifndef DWARFELEPHANTINITIALIZERBSYSTEMACTION_H
#define DWARFELEPHANTINITIALIZERBSYSTEMACTION_H

//libMesh includes
#include "libmesh/equation_systems.h"

#include "Action.h"
#include "MooseMesh.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"
#include "DisplacedProblem.h"

#include "DwarfElephantRBClasses.h"

// Forward Declarations
namespace libMesh
{
  class EquationSystems;
}

class MooseMesh;
class DisplacedProblem;

class DwarfElephantInitializeRBSystemAction : public Action
{
public:
  DwarfElephantInitializeRBSystemAction(InputParameters params);

  virtual void act() override;

  void initializeParameters();
  void initializeRBSystem();

protected:
  bool _use_displaced;
  bool _skip_matrix_assembly_in_rb_system;
  bool _skip_vector_assembly_in_rb_system;
  bool _offline_stage;

  std::string _parameters_filename;
  std::string _system_name;

  EquationSystems * _es_ptr;
  TransientNonlinearImplicitSystem * _sys;
  MooseMesh * _mesh_ptr;
  DwarfElephantRBConstruction * _rb_con_ptr;
};

template<>
InputParameters validParams<DwarfElephantInitializeRBSystemAction>();

#endif //DWARFELEPHANTINITIALIZERBSYSTEMACTION_H
