/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantConduction.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantConduction>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("The class implements a thermal conduction \
                              problem.");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantConduction::DwarfElephantConduction(const InputParameters & parameters) :
  Diffusion(parameters),
  // gets the thermal conductivity directly from the corresponding material file
  _lambda(getMaterialProperty<Real>("thermal_conductivity"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantConduction::computeQpResidual()
{
  return _lambda[_qp] * Diffusion::computeQpResidual();
}

Real
DwarfElephantConduction::computeQpJacobian()
{
  return _lambda[_qp] * Diffusion::computeQpJacobian();
}
