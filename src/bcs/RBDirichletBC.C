#include "RBDirichletBC.h"

template<>
InputParameters validParams<RBDirichletBC>()
{
  InputParameters p = validParams<RBNodalBC>();
  p.addRequiredParam<Real>("value", "Value of the BC");

  return p;
}


RBDirichletBC::RBDirichletBC(const InputParameters & parameters) :
  RBNodalBC(parameters),
  _value(getParam<Real>("value"))
{}

Real
RBDirichletBC::computeQpResidual()
{
  return _u[_qp] - _value;
}
