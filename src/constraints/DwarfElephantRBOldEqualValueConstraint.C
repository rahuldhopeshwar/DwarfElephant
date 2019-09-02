#include "DwarfElephantRBOldEqualValueConstraint.h"
#include "SubProblem.h"
#include "FEProblem.h"

registerMooseObject("MooseApp", DwarfElephantRBOldEqualValueConstraint);

template <>
InputParameters
validParams<DwarfElephantRBOldEqualValueConstraint>()
{
  InputParameters params = validParams<DwarfElephantRBMortarConstraint>();
  params.addClassDescription(
      "OldEqualValueConstraint enforces solution continuity between slave and "
      "master sides of a mortar interface using lagrange multipliers");
  return params;
}

DwarfElephantRBOldEqualValueConstraint::DwarfElephantRBOldEqualValueConstraint(const InputParameters & parameters)
  : DwarfElephantRBMortarConstraint(parameters)
{
}

Real
DwarfElephantRBOldEqualValueConstraint::computeQpResidual(Moose::MortarType mortar_type)
{
  switch (mortar_type)
  {
    case Moose::MortarType::Slave:
      return -_lambda[_qp] * _test_slave[_i][_qp];
    case Moose::MortarType::Master:
      return _lambda[_qp] * _test_master[_i][_qp];
    case Moose::MortarType::Lower:
      return (_u_master[_qp] - _u_slave[_qp]) * _test[_i][_qp];
    default:
      return 0;
  }
}

Real
DwarfElephantRBOldEqualValueConstraint::computeQpJacobian(Moose::ConstraintJacobianType jacobian_type,
                                           unsigned int jvar)
{
  typedef Moose::ConstraintJacobianType JType;

  switch (jacobian_type)
  {
    case JType::SlaveLower:
      if (jvar == _var->number())
        return -(*_phi)[_j][_qp] * _test_slave[_i][_qp];
      break;

    case JType::MasterLower:
      if (jvar == _var->number())
        return (*_phi)[_j][_qp] * _test_master[_i][_qp];
      break;

    case JType::LowerSlave:
      if (jvar == _slave_var.number())
        return -(*_phi)[_j][_qp] * _test[_i][_qp];
      break;

    case JType::LowerMaster:
      if (jvar == _master_var.number())
        return (*_phi)[_j][_qp] * _test[_i][_qp];
      break;

    default:
      return 0;
  }

  return 0;
}
