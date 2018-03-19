/**
 * This Kernel is implements a darcy flow using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEDarcyNormalized.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEDarcyNormalized>()
{
  InputParameters params = validParams<DwarfElephantFEDarcy>();

  params.addClassDescription("The class implements a darcy flow problem.");
  params.addRequiredParam<Real>("norm_value", "Defines the normalization value.");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEDarcyNormalized::DwarfElephantFEDarcyNormalized(const InputParameters & parameters) :
  DwarfElephantFEDarcy(parameters),
  _norm_value(getParam<Real>("norm_value"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEDarcyNormalized::computeQpResidual()
{
  return 1./_norm_value * DwarfElephantFEDarcy::computeQpResidual();
}

Real
DwarfElephantFEDarcyNormalized::computeQpJacobian()
{
  return 1./_norm_value * DwarfElephantFEDarcy::computeQpJacobian();
}
