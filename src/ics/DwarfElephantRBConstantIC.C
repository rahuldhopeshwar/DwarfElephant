#include "libmesh/point.h"

#include "DwarfElephantRBConstantIC.h"

template <>
InputParameters
validParams<DwarfElephantRBConstantIC>()
{
  InputParameters params = validParams<DwarfElephantRBInitialCondition>();
  params.addRequiredParam<Real>("value", "The value to be set in IC");

  params.addClassDescription("Sets a constant field value."
                             "This class is only required if you are using parameter-dependent "
                             "initial conditions.");
  return params;
}

DwarfElephantRBConstantIC::DwarfElephantRBConstantIC(const InputParameters & parameters)
  : DwarfElephantRBInitialCondition(parameters), _value(getParam<Real>("value"))
{
}

Real
DwarfElephantRBConstantIC::value(const Point & /*p*/)
{
  return _value;
}
