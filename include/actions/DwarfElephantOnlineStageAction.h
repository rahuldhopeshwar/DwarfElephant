#ifndef DWARFELEPHANTONLINESTAGEACTION_H
#define DWARFELEPHANTONLINESTAGEACTION_H

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

class DwarfElephantOnlineStageAction : public Action
{
public:
  DwarfElephantOnlineStageAction(InputParameters params);

  virtual void act() override;

};

template<>
InputParameters validParams<DwarfElephantOnlineStageAction>();

#endif //DWARFELEPHANTONLINESTAGEACTION_H
