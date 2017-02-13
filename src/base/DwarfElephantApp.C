/// MOOSE includes
#include "Moose.h"
#include "AppFactory.h"
//#include "ActionFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

/// MOOSE includes (DwarfElephant package)
#include "DwarfElephantApp.h"

// Kernels
#include "Conduction.h"
#include "RBKernel.h"
#include "RBDiffusion.h"

// AuxKernels
#include "KernelOutputAux.h"

// Materials
#include "SandStone.h"
#include "Shale.h"

//Outputs
#include "RBOutput.h"
#include "KernelOutput.h"

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
  // Kernels
  registerKernel(Conduction);
  registerKernel(RBKernel);
  registerKernel(RBDiffusion);

  // AuxKernels
  registerAux(KernelOutputAux);

  // Materials
  registerMaterial(SandStone);
  registerMaterial(Shale);

  // Outputs
  registerOutput(RBOutput);
  registerOutput(KernelOutput);
}

// External entry point for dynamic syntax association
extern "C" void DwarfElephantApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DwarfElephantApp::associateSyntax(syntax, action_factory); }
void
DwarfElephantApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
//    registerAction(KernelOutputAction, "add_postprocessor");
//    syntax.registerActionSyntax("KernelOutputAction", "KernelOutput");
}
