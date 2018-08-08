#include "DwarfElephantLiftingFunction.h"

template <>
InputParameters
validParams<DwarfElephantLiftingFunction>()
{
  InputParameters params = validParams<Function>();
  params.addClassDescription("This class implements the lifitng function for variable z-values.");
  params.addRequiredParam<FunctionName>("boundary_function", "The function that describes the value distribution of the boundary.");
  return params;
}

DwarfElephantLiftingFunction::DwarfElephantLiftingFunction(const InputParameters & parameters)
  : Function(parameters),
    FunctionInterface(this),
    _boundary_function(getFunctionByName("boundary_function"))
{
}

void
DwarfElephantLiftingFunction::initialSetup()
{
  _bottom_nodes = &_sc_fe_problem.mesh().getNodeList(5);

  // for(unsigned int i=_bottom_nodes[0];  i< _bottom_nodes->size(); i++)
  // {
  //   _console << _bottom_nodes[i] << std::endl;
  // }
}

Real
DwarfElephantLiftingFunction::value(Real /*t*/, const Point & /*p*/)
{

  Real value = 0.0;
  return value;
}
