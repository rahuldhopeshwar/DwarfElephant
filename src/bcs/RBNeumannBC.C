#include "RBNeumannBC.h"

template<>
InputParameters validParams<RBNeumannBC>()
{
  InputParameters params = validParams<RBIntegratedBC>();
  params.addParam<Real>("value", 0.0, "The value of the gradient on the boundary.");
  return params;
}

RBNeumannBC::RBNeumannBC(const InputParameters & parameters) :
  RBIntegratedBC(parameters),
  _value(getParam<Real>("value"))
{}

Real
RBNeumannBC::computeQpResidual()
{
  return -_test[_i][_qp]*_value;
}


