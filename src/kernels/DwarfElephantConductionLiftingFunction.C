/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 ///---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantConductionLiftingFunction.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantConductionLiftingFunction>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("The class implements a thermal conduction \
                              problem.");
  params.addRequiredParam<FunctionName>("lifting_function", "Name of the lifting function two account for the inhomogeneous Dirichlet boundary conditions.");		      
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantConductionLiftingFunction::DwarfElephantConductionLiftingFunction(const InputParameters & parameters) :
  Diffusion(parameters),
  // gets the thermal conductivity directly from the corresponding material file
  _lambda(getMaterialProperty<Real>("thermal_conductivity")),
  _lifting_function(&getFunction("lifting_function"))
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation

Real
DwarfElephantConductionLiftingFunction::computeQpResidual()
{
  return _lambda[_qp] *   (_grad_test[_i][_qp]*(_grad_u[_qp]-_lifting_function->gradient(_fe_problem.time(),_qp)));
}

Real
DwarfElephantConductionLiftingFunction::computeQpJacobian()
{
   return  _lambda[_qp]*(_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
}
