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
#include "RBNodalBC.h"
#include "RBDirichletBC.h"
#include "RBPresetNodalBC.h"
#include "RBPresetBC.h"
#include "RBIntegratedBC.h"
#include "RBNeumannBC.h"

// Kernels
#include "Conduction.h"
#include "DwarfElephantConductionLiftingFunction.h"
#include "Darcy.h"
#include "RBKernel.h"
#include "RBDiffusion.h"
#include "DwarfElephantRBDarcy.h"
#include "RBDiffusionLiftingFunction.h"

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
#include "DwarfElephantRBTransient.h"

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
  registerBoundaryCondition(RBNodalBC);
  registerBoundaryCondition(RBDirichletBC);
  registerBoundaryCondition(RBPresetNodalBC);
  registerBoundaryCondition(RBPresetBC);
  registerBoundaryCondition(RBIntegratedBC);
  registerBoundaryCondition(RBNeumannBC);

  // Kernels
  registerKernel(Conduction);
  registerKernel(DwarfElephantConductionLiftingFunction);
  registerKernel(Darcy);
  registerKernel(RBKernel);
  registerKernel(RBDiffusion);
  registerKernel(DwarfElephantRBDarcy);
  registerKernel(RBDiffusionLiftingFunction);

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
  registerExecutioner(DwarfElephantRBTransient);

}

// External entry point for dynamic syntax association
extern "C" void DwarfElephantApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DwarfElephantApp::associateSyntax(syntax, action_factory); }
void
DwarfElephantApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
