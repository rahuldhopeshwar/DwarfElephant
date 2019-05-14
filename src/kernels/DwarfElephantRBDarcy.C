/**
 * This Kernel is implements a darcy flow using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBDarcy.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBDarcy);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBDarcy>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();

  params.addClassDescription("The class implements a darcy flow problem.");
  params.addRequiredParam<Real>("permeability", "Defines the rock permeability.");
  params.addParam<Real>("fluid_viscosity", 1.0, "Defines the fluid viscosity.");
  params.addParam<Real>("norm_value_permeability", 1.0, "Defines the normalization value for the permeability.");
  params.addParam<Real>("norm_value_viscosity", 1.0, "Defines the normalization value of the fluid viscosity.");
  params.addParam<bool>("gravity_term", false, "Enables the gravity term.");
  params.addParam<Real>("fluid_density", 0.0, "Defines the fluid density.");
  params.addParam<RealVectorValue>("gravity", RealVectorValue(0.0, 0.0, 0.0),"Defines the gravity.");

  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBDarcy::DwarfElephantRBDarcy(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _permeability(getParam<Real>("permeability")),
  _norm_value_perm(getParam<Real>("norm_value_permeability")),
  _viscosity(getParam<Real>("fluid_viscosity")),
  _norm_value_visc(getParam<Real>("norm_value_viscosity")),
  _gravity_term(getParam<bool>("gravity_term")),
  _fluid_density(getParam<Real>("fluid_density")),
  _gravity(getParam<RealVectorValue>("gravity"))
{
  if((_gravity_term && _fluid_density == 0) || (_gravity_term && _gravity == RealVectorValue(0.0, 0.0, 0.0)))
    mooseError("When you use the gravity term you have to define a fluid density and a gravity");
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantRBDarcy::computeQpResidual()
{
  if(!_gravity_term)
    return (_permeability/_norm_value_perm) * (_norm_value_visc/_viscosity) * _grad_u[_qp] * _grad_test[_i][_qp];
  else
    return (_permeability/_norm_value_perm) * (_norm_value_visc/_viscosity) * (_grad_u[_qp]-(_fluid_density*_gravity)) * _grad_test[_i][_qp];
}

Real
DwarfElephantRBDarcy::computeQpJacobian()
{
  return (_permeability/_norm_value_perm) * (_norm_value_visc/_viscosity) * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
