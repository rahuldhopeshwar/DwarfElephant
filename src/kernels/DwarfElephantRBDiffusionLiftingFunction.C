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

///---------------------------------INCLUDES--------------------------------
//MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBDiffusionLiftingFunction.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBDiffusionLiftingFunction>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  params.addClassDescription("Implements a Diffusion problem using \
                             the RBKernel.");
  params.addRequiredParam<FunctionName>("lifting_function", "Name of the lifting function two account for the inhomogeneous Dirichlet boundary conditions.");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBDiffusionLiftingFunction::DwarfElephantRBDiffusionLiftingFunction(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _lifting_function(&getFunction("lifting_function"))
{
}

///----------------------------------PDEs-----------------------------------
Real
DwarfElephantRBDiffusionLiftingFunction::computeQpResidual()
{
  return    (_grad_test[_i][_qp]*(_grad_u[_qp]-_lifting_function->gradient(_fe_problem.time(),_qp)));
}

Real
DwarfElephantRBDiffusionLiftingFunction::computeQpJacobian()
{
  return  (_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
}

//Real
//DwarfElephantRBDiffusionLiftingFunction::computeQpMassMatrix()
//{
//  return _phi[_j][_qp] * _test[_i][_qp];
//}
