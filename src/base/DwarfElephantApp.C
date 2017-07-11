/// MOOSE includes
#include "Moose.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

/// MOOSE includes (DwarfElephant package)
#include "DwarfElephantApp.h"

// Base
#include "DwarfElephantRBProblem.h"

//BCs
#include "DwarfElephantRBNodalBC.h"
#include "DwarfElephantRBDirichletBC.h"
#include "DwarfElephantRBPresetNodalBC.h"
#include "DwarfElephantRBPresetBC.h"
#include "DwarfElephantRBIntegratedBC.h"
#include "DwarfElephantRBNeumannBC.h"

// Kernels
#include "DwarfElephantConduction.h"
#include "DwarfElephantConductionLiftingFunction.h"
#include "DwarfElephantDarcy.h"
#include "DwarfElephantRBKernel.h"
#include "DwarfElephantRBDiffusion.h"
#include "DwarfElephantRBDarcy.h"
#include "DwarfElephantRBDiffusionLiftingFunction.h"

// Materials
#include "SandStone.h"
#include "Shale.h"

// UserObjects
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"
#include "DwarfElephantOfflineOnlineStageSteadyState.h"
#include "DwarfElephantOfflineOnlineStageTransient.h"

// Functions
#include "CacheBoundaries.h"

// Executioners
#include "DwarfElephantRBSteady.h"

template<>
InputParameters validParams<DwarfElephantApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

DwarfElephantApp::DwarfElephantApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  DwarfElephantApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  DwarfElephantApp::associateSyntax(_syntax, _action_factory);
}

DwarfElephantApp::~DwarfElephantApp()
{
}

// External entry point for dynamic application loading
extern "C" void DwarfElephantApp__registerApps() { DwarfElephantApp::registerApps(); }
void
DwarfElephantApp::registerApps()
{
  registerApp(DwarfElephantApp);
}

// External entry point for dynamic object registration
extern "C" void DwarfElephantApp__registerObjects(Factory & factory) { DwarfElephantApp::registerObjects(factory); }
void
DwarfElephantApp::registerObjects(Factory & factory)
{
  // Base
  registerProblem(DwarfElephantRBProblem);

  // BCs
  registerBoundaryCondition(DwarfElephantRBNodalBC);
  registerBoundaryCondition(DwarfElephantRBDirichletBC);
  registerBoundaryCondition(DwarfElephantRBPresetNodalBC);
  registerBoundaryCondition(DwarfElephantRBPresetBC);
  registerBoundaryCondition(DwarfElephantRBIntegratedBC);
  registerBoundaryCondition(DwarfElephantRBNeumannBC);

  // Kernels
  registerKernel(DwarfElephantConduction);
  registerKernel(DwarfElephantConductionLiftingFunction);
  registerKernel(DwarfElephantDarcy);
  registerKernel(DwarfElephantRBKernel);
  registerKernel(DwarfElephantRBDiffusion);
  registerKernel(DwarfElephantRBDarcy);
  registerKernel(DwarfElephantRBDiffusionLiftingFunction);

  // Materials
  registerMaterial(SandStone);
  registerMaterial(Shale);

  // UserObjects
  registerUserObject(DwarfElephantInitializeRBSystemSteadyState);
  registerUserObject(DwarfElephantInitializeRBSystemTransient);
  registerUserObject(DwarfElephantOfflineOnlineStageSteadyState);
  registerUserObject(DwarfElephantOfflineOnlineStageTransient);

  // Functions
  registerFunction(CacheBoundaries);

  // Executioners
  registerExecutioner(DwarfElephantRBSteady);

}

// External entry point for dynamic syntax association
extern "C" void DwarfElephantApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DwarfElephantApp::associateSyntax(syntax, action_factory); }
void
DwarfElephantApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
