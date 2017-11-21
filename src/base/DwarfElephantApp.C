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
#include "DwarfElephantRBFunctionDirichletBC.h"
#include "DwarfElephantRBPenaltyDirichletBC.h"
#include "DwarfElephantRBPresetNodalBC.h"
#include "DwarfElephantRBPresetBC.h"
#include "DwarfElephantRBFunctionPresetBC.h"
#include "DwarfElephantRBIntegratedBC.h"
#include "DwarfElephantRBNeumannBC.h"
#include "DwarfElephantRBNeumannBCND.h"
#include "DwarfElephantRBFunctionNeumannBC.h"

// Kernels
#include "DwarfElephantFEThermalConduction.h"
#include "DwarfElephantFEConductionLiftingFunction.h"
#include "DwarfElephantFEElectricalConduction.h"
#include "DwarfElephantFEDarcy.h"
//#include "DwarfElephantRBKernel.h"
#include "DwarfElephantRBDiffusion.h"
#include "DwarfElephantRBDiffusionND.h"
#include "DwarfElephantRBDiffusionLiftingFunction.h"
#include "DwarfElephantRBTimeDerivative.h"
#include "ExtractQpPointsKernel.h"

// Materials
#include "DwarfElephantSandStone.h"
#include "DwarfElephantShale.h"
#include "DwarfElephantFault.h"
#include "ExtractQpPoints.h"

// UserObjects
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"
#include "DwarfElephantOfflineOnlineStageSteadyState.h"
#include "DwarfElephantOfflineOnlineStageTransient.h"
#include "DwarfElephantERTPreCalculations.h"

// Functions
#include "CacheBoundaries.h"

// Executioners
#include "DwarfElephantRBExecutioner.h"

// Outputs
#include "DwarfElephantDakotaOutput.h"
#include "DwarfElephantRBOutput.h"

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
  registerBoundaryCondition(DwarfElephantRBFunctionDirichletBC);
  registerBoundaryCondition(DwarfElephantRBPenaltyDirichletBC);
  registerBoundaryCondition(DwarfElephantRBPresetNodalBC);
  registerBoundaryCondition(DwarfElephantRBPresetBC);
  registerBoundaryCondition(DwarfElephantRBFunctionPresetBC);
  registerBoundaryCondition(DwarfElephantRBIntegratedBC);
  registerBoundaryCondition(DwarfElephantRBNeumannBC);
  registerBoundaryCondition(DwarfElephantRBNeumannBCND);
  registerBoundaryCondition(DwarfElephantRBFunctionNeumannBC);

  // Kernels
  registerKernel(DwarfElephantFEThermalConduction);
  registerKernel(DwarfElephantFEConductionLiftingFunction);
  registerKernel(DwarfElephantFEElectricalConduction);
  registerKernel(DwarfElephantFEDarcy);
//  registerKernel(DwarfElephantRBKernel);
  registerKernel(DwarfElephantRBDiffusion);
  registerKernel(DwarfElephantRBDiffusionND);
  registerKernel(DwarfElephantRBDiffusionLiftingFunction);
  registerKernel(DwarfElephantRBTimeDerivative);
  registerKernel(ExtractQpPointsKernel);

  // Materials
  registerMaterial(DwarfElephantSandStone);
  registerMaterial(DwarfElephantShale);
  registerMaterial(DwarfElephantFault);
  registerMaterial(ExtractQpPoints);

  // UserObjects
  registerUserObject(DwarfElephantInitializeRBSystemSteadyState);
  registerUserObject(DwarfElephantInitializeRBSystemTransient);
  registerUserObject(DwarfElephantOfflineOnlineStageSteadyState);
  registerUserObject(DwarfElephantOfflineOnlineStageTransient);
  registerUserObject(DwarfElephantERTPreCalculations);

  // Functions
  registerFunction(CacheBoundaries);

  // Executioners
  registerExecutioner(DwarfElephantRBExecutioner);

  // Outputs
  registerOutput(DwarfElephantDakotaOutput);
  registerOutput(DwarfElephantRBOutput);

}

// External entry point for dynamic syntax association
extern "C" void DwarfElephantApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DwarfElephantApp::associateSyntax(syntax, action_factory); }
void
DwarfElephantApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
