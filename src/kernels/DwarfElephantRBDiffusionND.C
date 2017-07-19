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
#include "DwarfElephantRBDiffusionND.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBDiffusionND>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  params.addClassDescription("Implements a Diffusion problem using \
                             the RBKernel.");
  params.addRequiredParam<Real>("u_ref", "Reference value for non dimensionalizing the problem");
  params.addRequiredParam<Real>("l_ref", "Reference length");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBDiffusionND::DwarfElephantRBDiffusionND(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _u_ref(getParam<Real>("u_ref")),
  _l_ref(getParam<Real>("l_ref"))
{
}

///----------------------------------PDEs-----------------------------------
Real
DwarfElephantRBDiffusionND::computeQpResidual()
{
  return (1/_u_ref) * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
DwarfElephantRBDiffusionND::computeQpJacobian()
{
  return (1/_u_ref) * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

Real
DwarfElephantRBDiffusionND::computeQpMassMatrix()
{
  return _phi[_j][_qp] * _test[_i][_qp];
}
