/**
 * This Kernel is implements a darcy flow using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "Darcy.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<Darcy>()
{
  InputParameters params = validParams<Diffusion>();

  params.addClassDescription("The class implements a darcy flow problem.");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
Darcy::Darcy(const InputParameters & parameters) :
  Diffusion(parameters),
  _permeability(getMaterialProperty<Real>("permeability")),
  _dynamic_viscosity(getMaterialProperty<Real>("dynamic_viscosity")),
  _fluid_density(getMaterialProperty<Real>("fluid_density")),
  _gravity(getMaterialProperty<RealVectorValue>("gravity"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
Darcy::computeQpResidual()
{
  return (_permeability[_qp]/_dynamic_viscosity[_qp]) * Diffusion::computeQpResidual() - ((_fluid_density[_qp] * _gravity[_qp]) * _grad_test[_i][_qp]);
}

Real
Darcy::computeQpJacobian()
{
  return (_permeability[_qp]/_dynamic_viscosity[_qp]) * Diffusion::computeQpJacobian();
}
