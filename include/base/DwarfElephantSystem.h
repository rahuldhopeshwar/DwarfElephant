#ifndef DWARFELEPHANTSYSTEM_H
#define DWARFELEPHANTSYSTEM_H

// libMesh includes (RB package)
#include "libmesh/rb_construction.h"

// MOOSE includes
#include "FEProblemBase.h"
#include "TimeIntegrator.h"
#include "NonlinearSystem.h"
#include "DwarfElephantRBAssembly.h"

// MOOSE includes (DwarfElephant package)
//#include "RBStructuresP1Theta3ThetaEqualMuSteadyState.h"
//#include "RBStructuresP1Theta5ThetaEqualMuSteadyState.h"

// Forward Declarations
namespace libMesh
{
  class RBConstruction;
}

class FEProblemBase;

class DwarfElephantSystem : public NonlinearSystem

{
public:
  DwarfElephantSystem(FEProblemBase & problem, const std::string & name);

  ~DwarfElephantSystem();

  // Initialize data structure
  virtual void 	solve () override;

  unsigned int u_var;

  DwarfElephantRBAssembly * _rb_assembly;
};

#endif /* DWARFELEPHANTSYSTEM_H */
