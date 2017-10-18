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
  params.addRequiredParam<Real>("thermal_conductivity", "Defines the value of the thermal conductivity");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEThermalConduction::DwarfElephantFEThermalConduction(const InputParameters & parameters) :
  Diffusion(parameters),
  // gets the thermal conductivity directly from the corresponding material file for future use in case of varying parameters
//  _lambda(getMaterialProperty<Real>("thermal_conductivity"))
  _lambda(getParam<Real>("thermal_conductivity"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEThermalConduction::computeQpResidual()
{
//  return _lambda[_qp] * Diffusion::computeQpResidual(); // in case of getting lambda from the material property class
  return _lambda * Diffusion::computeQpResidual();
}

Real
DwarfElephantFEThermalConduction::computeQpJacobian()
{
//  return _lambda[_qp] * Diffusion::computeQpJacobian(); // in case of getting lambda from the material property class
  return _lambda * Diffusion::computeQpJacobian();
}
