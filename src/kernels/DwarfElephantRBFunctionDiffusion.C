/**
 * This Kernel implements the Diffusion problem by using the RB method.
 * It is important to note that every PDE that will be used within the RB
 * method has to inherit from RBKernel and not from Kernel.
 * Furthermore, one should recall that the stiffness matrix and the load
 * vector are constructed out of the parameter independent part of the
 * PDE. Therefore, the functions computeQpResidual() and computeQpJacobian()
 * should only contain this parameter independent part. Consequently, the
 * RBDiffusion Kernel is used for a Conduction problem.
 */

//---------------------------------INCLUDES--------------------------------
//MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBFunctionDiffusion.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBFunctionDiffusion);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBFunctionDiffusion>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  params.addClassDescription("Implements a Diffusion problem using \
                             the RBKernel.");
  params.addRequiredParam<FunctionName>("function", "Name of the function that describes the parameter dependence.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBFunctionDiffusion::DwarfElephantRBFunctionDiffusion(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _func(getFunction("function"))
{
}

//----------------------------------PDEs-----------------------------------
Real
DwarfElephantRBFunctionDiffusion::computeQpResidual()
{
  Real dependency = _func.value(_fe_problem.time(),_q_point[_qp]);
  return -dependency*(_grad_u[_qp] * _grad_test[_i][_qp]);
}

Real
DwarfElephantRBFunctionDiffusion::computeQpJacobian()
{
  Real dependency = _func.value(_fe_problem.time(),_q_point[_qp]);
  return dependency*_grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

Real
DwarfElephantRBFunctionDiffusion::computeQpOutput()
{
  return 1;
}
