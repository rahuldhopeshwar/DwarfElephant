#include "ConductionThreeLayerApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

// Actions
#include "RBSimpleConstruction.h"

// Kernels
#include "Conduction.h"
#include "RBKernel.h"

// Materials
#include "SandStone.h"
#include "Shale.h"
#include "RBAssembly.h"

// UserObjects
#include "GeneralizedPlaneStrainUserObject.h"

//Outputs
#include "RBOutput.h"

template<>
InputParameters validParams<ConductionThreeLayerApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

ConductionThreeLayerApp::ConductionThreeLayerApp(InputParameters parameters) :
    MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ConductionThreeLayerApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ConductionThreeLayerApp::associateSyntax(_syntax, _action_factory);
}

ConductionThreeLayerApp::~ConductionThreeLayerApp()
{
}

void
ConductionThreeLayerApp::registerObjects(Factory & factory)
{
  // Register any custom objects you have built on the MOOSE Framework
  // Kernels
  registerKernel(Conduction);
  registerKernel(RBKernel);

  // Materials
  registerMaterial(SandStone);
  registerMaterial(Shale);
  registerMaterial(RBAssembly);

  // UserObjects
  registerUserObject(GeneralizedPlaneStrainUserObject);

  // Outputs
  registerOutput(RBOutput);
}

void
ConductionThreeLayerApp::registerApps()
{
  registerApp(ConductionThreeLayerApp);
}

void
ConductionThreeLayerApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  registerAction(RBSimpleConstruction, "add_user_object");
  syntax.registerActionSyntax("RBSimpleConstruction", "RBSimpleConstruction");
}

