/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEThermalConduction.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEThermalConduction>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("The class implements a thermal conduction \
                              problem.");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEThermalConduction::DwarfElephantFEThermalConduction(const InputParameters & parameters) :
  Diffusion(parameters),
  // gets the thermal conductivity directly from the corresponding material file
  _lambda(getMaterialProperty<Real>("thermal_conductivity"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEThermalConduction::computeQpResidual()
{
  return _lambda[_qp] * Diffusion::computeQpResidual();
}

Real
DwarfElephantFEThermalConduction::computeQpJacobian()
{
  return _lambda[_qp] * Diffusion::computeQpJacobian();
}
