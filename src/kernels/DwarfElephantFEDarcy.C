/**
 * This Kernel is implements a darcy flow using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEDarcy.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEDarcy>()
{
  InputParameters params = validParams<Diffusion>();

  params.addClassDescription("The class implements a darcy flow problem.");
  params.addRequiredParam<Real>("permeability", "Defines the rock permeability in relationship to a certain reference permeability.");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEDarcy::DwarfElephantFEDarcy(const InputParameters & parameters) :
  Diffusion(parameters),
  _permeability(getParam<Real>("permeability"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEDarcy::computeQpResidual()
{
  return _permeability * Diffusion::computeQpResidual();
}

Real
DwarfElephantFEDarcy::computeQpJacobian()
{
  return _permeability * Diffusion::computeQpJacobian();
}
