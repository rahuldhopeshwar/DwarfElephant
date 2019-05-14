/**
 * This Kernel is implements a constant radiogenic heat production using the
 * Reduced Basis solution.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBVariableTimeDerivative.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBVariableTimeDerivative);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBVariableTimeDerivative>()
{
  InputParameters params = validParams<DwarfElephantRBTimeDerivative>();
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBVariableTimeDerivative::DwarfElephantRBVariableTimeDerivative(const InputParameters & parameters) :
  DwarfElephantRBTimeDerivative(parameters)
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantRBVariableTimeDerivative::computeQpResidual()
{
  return  _u_dot[_qp] * _test[_i][_qp];
}

Real
DwarfElephantRBVariableTimeDerivative::computeQpJacobian()
{
  return _du_dot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}
