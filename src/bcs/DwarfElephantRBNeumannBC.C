#include "DwarfElephantRBNeumannBC.h"

template<>
InputParameters validParams<DwarfElephantRBNeumannBC>()
{
  InputParameters params = validParams<DwarfElephantRBIntegratedBC>();
  params.addParam<Real>("value", 0.0, "The value of the gradient on the boundary.");
  return params;
}

DwarfElephantRBNeumannBC::DwarfElephantRBNeumannBC(const InputParameters & parameters) :
  DwarfElephantRBIntegratedBC(parameters),
  _value(getParam<Real>("value"))
{}

Real
DwarfElephantRBNeumannBC::computeQpResidual()
{
  return _test[_i][_qp]*_value;
}
