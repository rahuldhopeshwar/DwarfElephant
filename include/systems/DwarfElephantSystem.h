/**
 * This System class is required to use the RB system instead of the FE
 * system. Furthermore, it prevents an additional solve.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTSYSTEM_H
#define DWARFELEPHANTSYSTEM_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes (RB package)
#include "libmesh/rb_construction.h"

// MOOSE includes
#include "FEProblemBase.h"
#include "TimeIntegrator.h"
#include "NonlinearSystem.h"
#include "KernelBase.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBClassesSteadyState.h"
#include "DwarfElephantRBKernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class RBConstruction;
}

class FEProblemBase;
class KernelBase;
class DwarfElephantRBConstructionSteadyState;

///This System class is required to use the RB system instead of the FE system. Furthermore, it prevents an additional solve.
class DwarfElephantSystem : public NonlinearSystem

{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantSystem(FEProblemBase & problem, const std::string & name);

  ~DwarfElephantSystem();

  /* Methods */
  // Initialize data structure
  // virtual void initialSetup() override;
  virtual void 	solve () override;
};

#endif /* DWARFELEPHANTSYSTEM_H */
