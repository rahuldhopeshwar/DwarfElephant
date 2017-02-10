/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "Conduction.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<Conduction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The class implements a thermal conduction \
                              problem.");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
Conduction::Conduction(const InputParameters & parameters) :
  Kernel(parameters),
  // gets the thermal conductivity directly from the corresponding material file
  _lambda(getMaterialProperty<Real>("conductivity"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
Conduction::computeQpResidual()
{
  return _lambda[_qp]*_grad_u[_qp] * _grad_test[_i][_qp];
}

Real
Conduction::computeQpJacobian()
{
  return _lambda[_qp]*_grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
