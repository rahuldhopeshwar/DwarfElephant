#ifndef DWARFELEPHANTSYSTEM_H
#define DWARFELEPHANTSYSTEM_H

// libMesh includes (RB package)
#include "libmesh/rb_construction.h"

// MOOSE includes
#include "FEProblemBase.h"
#include "TimeIntegrator.h"
#include "NonlinearSystem.h"

// MOOSE includes (DwarfElephant package)
//#include "RBStructuresP1Theta3ThetaEqualMuSteadyState.h"
//#include "RBStructuresP1Theta5ThetaEqualMuSteadyState.h"
#include "DwarfElephantRBClassesSteadyState.h"

// Forward Declarations
namespace libMesh
{
  class RBConstruction;
}

class FEProblemBase;
class DwarfElephantRBConstructionSteadyState;

class DwarfElephantSystem : public NonlinearSystem

{
public:
  DwarfElephantSystem(FEProblemBase & problem, const std::string & name);

  ~DwarfElephantSystem();

  // Initialize data structure
  virtual void 	solve () override;
};

#endif /* DWARFELEPHANTSYSTEM_H */
