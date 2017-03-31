#include "RBDirichletBC.h"

template<>
InputParameters validParams<RBDirichletBC>()
{
  InputParameters p = validParams<RBNodalBC>();
  p.addRequiredParam<Real>("value", "Value of the BC");
  p.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");

  return p;
}


RBDirichletBC::RBDirichletBC(const InputParameters & parameters) :
  RBNodalBC(parameters),
  _value(getParam<Real>("value")),
  _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject"))
{}

Real
RBDirichletBC::computeQpResidual()
{
  return _u[_qp] - _value;

//  _initialize_rb_system._solution->set(_current_node, -value);
}
