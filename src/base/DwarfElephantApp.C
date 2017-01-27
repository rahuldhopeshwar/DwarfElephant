#include "DwarfElephantApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
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

//Outputs
#include "RBOutput.h"

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
//  // Register any custom objects you have built on the MOOSE Framework
//  // Kernels
  registerKernel(Conduction);
  registerKernel(RBKernel);

  // Materials
  registerMaterial(SandStone);
  registerMaterial(Shale);
  registerMaterial(RBAssembly);

  // Outputs
  registerOutput(RBOutput);
}

// External entry point for dynamic syntax association
extern "C" void DwarfElephantApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DwarfElephantApp::associateSyntax(syntax, action_factory); }
void
DwarfElephantApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  registerAction(RBSimpleConstruction, "add_user_object");
  syntax.registerActionSyntax("RBSimpleConstruction", "RBSimpleConstruction");
}
