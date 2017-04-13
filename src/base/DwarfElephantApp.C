/// MOOSE includes
#include "Moose.h"
#include "AppFactory.h"
//#include "ActionFactory.h"
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

// Kernels
#include "Conduction.h"
#include "RBKernel.h"
#include "RBDiffusion.h"

// Materials
#include "SandStone.h"
#include "Shale.h"

// UserObjects
#include "DwarfElephantInitializeRBSystem.h"
#include "DwarfElephantOfflineOnlineStage.h"

// Functions
#include "CacheBoundaries.h"

// Executioners
#include "DwarfElephantExecutioner.h"


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

  // Kernels
  registerKernel(Conduction);
  registerKernel(RBKernel);
  registerKernel(RBDiffusion);

  // Materials
  registerMaterial(SandStone);
  registerMaterial(Shale);

  // UserObjects
  registerUserObject(DwarfElephantInitializeRBSystem);
  registerUserObject(DwarfElephantOfflineOnlineStage);

  // Functions
  registerFunction(CacheBoundaries);

  // Executioners
  registerExecutioner(DwarfElephantExecutioner);

}

// External entry point for dynamic syntax association
extern "C" void DwarfElephantApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DwarfElephantApp::associateSyntax(syntax, action_factory); }
void
DwarfElephantApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
//    registerAction(DwarfElephantOnlineStageAction, "add_kernel");
//    syntax.registerActionSyntax("DwarfElephantOnlineStageAction", "OnlineStage");
}
