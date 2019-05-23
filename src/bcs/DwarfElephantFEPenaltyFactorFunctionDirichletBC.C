#include "DwarfElephantFEPenaltyFactorFunctionDirichletBC.h"
#include "Function.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFEPenaltyFactorFunctionDirichletBC);

template <>
InputParameters
validParams<DwarfElephantFEPenaltyFactorFunctionDirichletBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("value", "Penalty scalar");
  params.addRequiredParam<FunctionName>("function", "Forcing function");

  return params;
}

DwarfElephantFEPenaltyFactorFunctionDirichletBC::DwarfElephantFEPenaltyFactorFunctionDirichletBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _func(getFunction("function")), _value(getParam<Real>("value"))
{
}

Real
DwarfElephantFEPenaltyFactorFunctionDirichletBC::computeQpResidual()
{
  return _func.value(_t, _q_point[_qp]) * _test[_i][_qp] * (-_value + _u[_qp]);
}

Real
DwarfElephantFEPenaltyFactorFunctionDirichletBC::computeQpJacobian()
{
  return _func.value(_t, _q_point[_qp]) * _phi[_j][_qp] * _test[_i][_qp];
}
