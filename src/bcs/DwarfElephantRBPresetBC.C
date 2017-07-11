#include "DwarfElephantRBPresetBC.h"

template<>
InputParameters validParams<DwarfElephantRBPresetBC>()
{
  InputParameters p = validParams<DwarfElephantRBNodalBC>();
  p.addRequiredParam<Real>("value", "Value of the BC");
  return p;
}


DwarfElephantRBPresetBC::DwarfElephantRBPresetBC(const InputParameters & parameters) :
  DwarfElephantRBPresetNodalBC(parameters),
  _value(getParam<Real>("value"))
{

}

Real
DwarfElephantRBPresetBC::computeQpValue()
{
  return _value;
}

