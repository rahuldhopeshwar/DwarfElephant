/**
 * This Kernel is implements a darcy flow problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes. It was slightly modified from the original DwarfElephantFEDarcy
 * Kernel in order to allow an easier usability within the Blackbox Wrapper of
 * OpenDA.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEDarcyOpenDA.h"
#include "MooseMesh.h"

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEDarcyOpenDA>()
{
  InputParameters params = validParams<Kernel>();

  params.addClassDescription("The class implements a darcy flow problem.");
  params.addRequiredParam<std::vector<Real>>("permeability", "Defines the rock permeability in relationship to a certain reference permeability.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEDarcyOpenDA::DwarfElephantFEDarcyOpenDA(const InputParameters & parameters) :
  Kernel(parameters),
  _permeability(getParam<std::vector<Real>>("permeability")),
  _ID_first_block(*_fe_problem.mesh().meshSubdomains().begin())
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEDarcyOpenDA::computeQpResidual()
{
  return _permeability[_current_elem->subdomain_id() - _ID_first_block] * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
DwarfElephantFEDarcyOpenDA::computeQpJacobian()
{
  return _permeability[_current_elem->subdomain_id() - _ID_first_block] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
