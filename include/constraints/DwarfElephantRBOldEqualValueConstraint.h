#pragma once

#include "DwarfElephantRBMortarConstraint.h"

class DwarfElephantRBOldEqualValueConstraint;

template <>
InputParameters validParams<DwarfElephantRBOldEqualValueConstraint>();

/**
 * Constrain the value of a variable to be the same on both sides of an
 * interface.
 */
class DwarfElephantRBOldEqualValueConstraint : public DwarfElephantRBMortarConstraint
{
public:
  DwarfElephantRBOldEqualValueConstraint(const InputParameters & parameters);

protected:
  Real computeQpResidual(Moose::MortarType mortar_type) final;
  Real computeQpJacobian(Moose::ConstraintJacobianType jacobian_type, unsigned int jvar) final;
};
