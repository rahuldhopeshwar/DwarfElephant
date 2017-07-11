#include "DwarfElephantRBDirichletBC.h"

template<>
InputParameters validParams<DwarfElephantRBDirichletBC>()
{
  InputParameters p = validParams<DwarfElephantRBNodalBC>();
  p.addRequiredParam<Real>("value", "Value of the BC");

  return p;
}


DwarfElephantRBDirichletBC::DwarfElephantRBDirichletBC(const InputParameters & parameters) :
  DwarfElephantRBNodalBC(parameters),
  _value(getParam<Real>("value"))
{}

Real
DwarfElephantRBDirichletBC::computeQpResidual()
{
  return _u[_qp] - _value;
}
