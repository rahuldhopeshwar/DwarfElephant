/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEThermalConductionNormalized.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEThermalConductionNormalized>()
{
  InputParameters params = validParams<DwarfElephantFEThermalConduction>();
  params.addClassDescription("The class implements a thermal conduction \
                              problem.");
  params.addRequiredParam<Real>("norm_value", "Defines the normalization value.");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEThermalConductionNormalized::DwarfElephantFEThermalConductionNormalized(const InputParameters & parameters) :
  DwarfElephantFEThermalConduction(parameters),
  _norm_value(getParam<Real>("norm_value"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEThermalConductionNormalized::computeQpResidual()
{
  return 1./_norm_value * DwarfElephantFEThermalConduction::computeQpResidual();
}

Real
DwarfElephantFEThermalConductionNormalized::computeQpJacobian()
{
  return 1./_norm_value * DwarfElephantFEThermalConduction::computeQpJacobian();
}
