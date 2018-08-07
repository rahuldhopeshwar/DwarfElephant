#include "DwarfElephantReverseLiftingFunction.h"

// MOOSE includes
#include "LineSegment.h"
#include "RayTracing.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<DwarfElephantReverseLiftingFunction>()
{
  InputParameters params = validParams<NodalVectorPostprocessor>();
  params.addRequiredParam<FunctionName>("lifting_function", "The lifting function that should be reversed.");
  params.addParam<std::string>("system", "rb0", "The name of the used system.");

  return params;
}

DwarfElephantReverseLiftingFunction::DwarfElephantReverseLiftingFunction(const InputParameters & parameters)
  : NodalVectorPostprocessor(parameters),
  _lifting_function(getFunction("lifting_function")),
  _nodal_solution_original(declareVector("nodal_solution_original")),
  _system(getParam<std::string>("system"))
{
}

void
DwarfElephantReverseLiftingFunction::initialize()
{
  if(!_unique_node_execute)
    mooseError("This VectorPostprocessor has to be execute with 'unique_node_execute' equal to 'true'.");

  // Initialize solution vector
  _nodal_solution = _fe_problem.es().get_system<NonlinearImplicitSystem>(_system).current_local_solution.get();

  _nodal_solution_original.clear();
}

void
DwarfElephantReverseLiftingFunction::execute()
{
  Point _point(_current_node->operator()(0), _current_node->operator()(1), _current_node->operator()(2));
  _nodal_solution_original.push_back(_nodal_solution->el(_current_node->id()));
  _nodal_solution->set(_current_node->id(),_nodal_solution->el(_current_node->id()) + _lifting_function.value(_t, _point));
}

void
DwarfElephantReverseLiftingFunction::finalize()
{
  _nodal_solution->close();
  *_fe_problem.es().get_system(_system).solution = *_nodal_solution;
}
