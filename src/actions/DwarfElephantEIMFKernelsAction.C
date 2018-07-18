//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DwarfElephantEIMFKernelsAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

//registerMooseAction("ExampleApp", DwarfElephantEIMFKernelsAction, "add_kernel");

template <>
InputParameters
validParams<DwarfElephantEIMFKernelsAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<NonlinearVariableName>("variable", "The name of the variable in the simulation");
  params.addRequiredParam<UserObjectName>("initialize_rb_system","The name of the DwarfElephantInitializeRBSystemSteadyState object");

  return params;
}

DwarfElephantEIMFKernelsAction::DwarfElephantEIMFKernelsAction(InputParameters params) : Action(params)
 {}

void
DwarfElephantEIMFKernelsAction::act()
{
  NonlinearVariableName variable =
      getParam<NonlinearVariableName>("variable");
  const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = _problem -> getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
  UserObjectName _initial_rb_system_name = getParam<UserObjectName>("initialize_rb_system");
  unsigned int _n_eim_basis_functions = _initialize_rb_system._eim_con_ptr -> get_rb_evaluation().get_n_basis_functions();

  // Do some error checking
  mooseAssert(variable.size() == 1, "Expected 1 variable, received " << variable.size());

  for (unsigned int _i_eim_basis_functions = 0; _i_eim_basis_functions < _n_eim_basis_functions; _i_eim_basis_functions++)
  {
    InputParameters params = _factory.getValidParams("DwarfElephantEIMFKernel");
    params.set<NonlinearVariableName>("variable") = variable;
	params.set<unsigned int>("eim_basis_function_index") = _i_eim_basis_functions;
	params.set<UserObjectName>("initialize_rb_system") = getParam<UserObjectName>("initialize_rb_system");//_initial_rb_system_name; 
    _problem->addKernel("DwarfElephantEIMFKernel", "EIM_FKernel_" + std::to_string(_i_eim_basis_functions), params);
  }
}
