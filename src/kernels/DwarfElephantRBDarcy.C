/**
 * This Kernel is implements a darcy flow using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBDarcy.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBDarcy>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();

  params.addClassDescription("The class implements a RB darcy flow problem.");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBDarcy::DwarfElephantRBDarcy(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _permeability(getMaterialProperty<Real>("permeability")),
  _dynamic_viscosity(getMaterialProperty<Real>("dynamic_viscosity")),
  _fluid_density(getMaterialProperty<Real>("fluid_density")),
  _gravity(getMaterialProperty<RealVectorValue>("gravity"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantRBDarcy::computeQpResidual()
{
  return (_permeability[_qp]/_dynamic_viscosity[_qp]) * ((_grad_u[_qp] * _grad_test[_i][_qp]) - ((_fluid_density[_qp] * _gravity[_qp]) * _grad_test[_i][_qp]));
}

Real
DwarfElephantRBDarcy::computeQpJacobian()
{
  return (_permeability[_qp]/_dynamic_viscosity[_qp]) * (_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
}
