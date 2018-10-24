#include "DwarfElephantRBNeumannBCND.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBNeumannBCND);

template<>
InputParameters validParams<DwarfElephantRBNeumannBCND>()
{
  InputParameters params = validParams<DwarfElephantRBIntegratedBC>();
  params.addParam<Real>("value", 0.0, "The value of the gradient on the boundary.");
  params.addRequiredParam<Real>("u_ref", "Reference value for non dimensionalizing the problem");
  params.addRequiredParam<Real>("l_ref", "Reference length");
  return params;
}

DwarfElephantRBNeumannBCND::DwarfElephantRBNeumannBCND(const InputParameters & parameters) :
  DwarfElephantRBIntegratedBC(parameters),
  _value(getParam<Real>("value")),
  _u_ref(getParam<Real>("u_ref")),
  _l_ref(getParam<Real>("l_ref"))
{}

Real
DwarfElephantRBNeumannBCND::computeQpResidual()
{
  return -(_l_ref/_u_ref)*_test[_i][_qp]*_value;
}
