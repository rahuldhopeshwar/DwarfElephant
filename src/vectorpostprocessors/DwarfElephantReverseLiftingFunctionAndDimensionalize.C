#include "DwarfElephantReverseLiftingFunctionAndDimensionalize.h"

// MOOSE includes
#include "LineSegment.h"
#include "RayTracing.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<DwarfElephantReverseLiftingFunctionAndDimensionalize>()
{
  InputParameters params = validParams<NodalVectorPostprocessor>();
  params.addRequiredParam<FunctionName>("lifting_function", "The lifting function that should be reversed.");
  params.addParam<std::string>("system", "rb0", "The name of the used system.");
  params.addParam<Real>("reference_value_variable", 1.0, "The reference value used for the nondimensionalization.");
  params.addParam<bool>("dimensionalize", true, "When set to true the nondimensionalizated values are backtransformed to their original values.");
  params.addParam<bool>("scale_and_add", true, "When true: new_variable = old_variable*ref+ref, when false: new_variable = old_variable*ref");
  params.addParam<bool>("reverse_lifting_function", true, "Determines whether a lifting function is reversed or not.");
  return params;
}

DwarfElephantReverseLiftingFunctionAndDimensionalize::DwarfElephantReverseLiftingFunctionAndDimensionalize(const InputParameters & parameters)
  : NodalVectorPostprocessor(parameters),
  _lifting_function(getFunction("lifting_function")),
  _nodal_solution_original(declareVector("nodal_solution_original")),
  _system(getParam<std::string>("system")),
  _reference_value_variable(getParam<Real>("reference_value_variable")),
  _dimensionalize(getParam<bool>("dimensionalize")),
  _scale_and_add(getParam<bool>("scale_and_add")),
  _reverse_lifting_function(getParam<bool>("reverse_lifting_function"))
{
}

void
DwarfElephantReverseLiftingFunctionAndDimensionalize::initialize()
{
  if(!_unique_node_execute)
    mooseError("This VectorPostprocessor has to be execute with 'unique_node_execute' equal to 'true'.");

  // Initialize solution vector
  _nodal_solution = _fe_problem.es().get_system<NonlinearImplicitSystem>(_system).current_local_solution.get();

  _nodal_solution_original.clear();
}

void
DwarfElephantReverseLiftingFunctionAndDimensionalize::execute()
{
  // store the original variable values in a VectorPostprocessor
  _nodal_solution_original.push_back(_nodal_solution->el(_current_node->id()));

  Real _value = _nodal_solution->el(_current_node->id());

  // reverse the lifitng function
  if(_reverse_lifting_function)
  {
    // Define a point for the lifting function
    Point _point(_current_node->operator()(0), _current_node->operator()(1), _current_node->operator()(2));
    // _value += _nodal_solution->el(_current_node->id()) + _lifting_function.value(_t, _point);
    _value += _lifting_function.value(_t, _point);
  }

  // dimensionalize the variable
  if(_dimensionalize)
  {
    _value *= _reference_value_variable;
    if(_scale_and_add)
      _value += _reference_value_variable;
  }

  // set the value
  _nodal_solution->set(_current_node->id(), _value);
}

void
DwarfElephantReverseLiftingFunctionAndDimensionalize::finalize()
{
  _nodal_solution->close();
  *_fe_problem.es().get_system(_system).solution = *_nodal_solution;
}
