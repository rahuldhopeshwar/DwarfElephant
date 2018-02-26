#include "DwarfElephantAllElementalVariableValues.h"

// MOOSE includes
#include "LineSegment.h"
#include "RayTracing.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<DwarfElephantAllElementalVariableValues>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addRequiredParam<VariableName>("variable", "The variable you want to be monitor.");

  return params;
}

DwarfElephantAllElementalVariableValues::DwarfElephantAllElementalVariableValues(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _elem_values(declareVector("elem_values")),
    _var_name(parameters.get<VariableName>("variable"))
{
}

void
DwarfElephantAllElementalVariableValues::initialize()
{
  _elem_values.clear();
}

void
DwarfElephantAllElementalVariableValues::execute()
{
  Real value = 0;
  unsigned int num_elems = _fe_problem.mesh().nElem();

  _elem_values.resize(num_elems);

  // Get Element Values
  for (unsigned int i = 0; i < num_elems; i++){
    Elem * _element = _fe_problem.mesh().getMesh().query_elem_ptr(i);
    _subproblem.prepare(_element, _tid);
    _subproblem.reinitElem(_element, _tid);
    MooseVariable & var = _subproblem.getVariable(_tid, _var_name);
    const VariableValue & u = var.sln();

    unsigned int n = u.size();
    for (unsigned int i = 0; i < n; i++)
      value += u[i];

    value /= n;
    _elem_values[i] = value;
  }
}

// void
// DwarfElephantAllElementalVariableValues::threadJoin(const UserObject & s){
//
// }
