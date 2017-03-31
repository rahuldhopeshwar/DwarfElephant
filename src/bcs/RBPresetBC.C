#include "RBPresetBC.h"

template<>
InputParameters validParams<RBPresetBC>()
{
  InputParameters p = validParams<RBNodalBC>();
  p.addRequiredParam<Real>("value", "Value of the BC");
  return p;
}


RBPresetBC::RBPresetBC(const InputParameters & parameters) :
  RBPresetNodalBC(parameters),
  _value(getParam<Real>("value"))
{

}

Real
RBPresetBC::computeQpValue()
{
  return _value;
}

