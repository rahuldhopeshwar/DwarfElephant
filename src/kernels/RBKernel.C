///---------------------------------INCLUDE---------------------------------
// MOOSE includes (own)
#include "RBKernel.h"

///-------------------------------------------------------------------------
template<>

///----------------------------INPUT PARAMETERS-----------------------------
InputParameters validParams<RBKernel>()
{
  // General
  InputParameters params = validParams<Kernel>();

  // For the parameters obtained out of the Input file
  params.addClassDescription("This class implements the Reduced Basis \
		                     Method (a Model Order Reduction technique).");
  params.addRequiredParam<unsigned int>("muTest", "");
  params.addParam<unsigned int>("test",20,"");
  return params;
}

///-------------------------------------------------------------------------
RBKernel::RBKernel(const InputParameters & parameters) :

  // General
  Kernel(parameters),
  _mu0(getMaterialProperty<RealVectorValue>("RBParameter0")),
  _muTest(getParam<unsigned int>("muTest")),
  _test(getParam<unsigned int>("test"))
{
}

///-----------------------------------PDEs----------------------------------
// THE PDE SHOULD BE IMPLEMENTED IN RBCONDUCTION KERNEL
// ONLY PARAMTER INDEPENDENT PART
Real
RBKernel::computeQpResidual()
{
  return _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
RBKernel::computeQpJacobian()
{
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

///-----------------------------------GET-----------------------------------
DenseMatrix<Number> &
RBKernel::getJacobian()
{
  return _local_ke;
}

DenseVector<Number> &
RBKernel::getResidual()
{
  return _local_re;
}
///--------------------------------------------------------------------------
//struct A0 : ElemAssembly
//{
//  RBKernel * _rb_ptr = libmesh_nullptr;
//  _rb_ptr->computeJacobian();
//}
////#############################################################################################################################################
//// In the case of a general RBKernel use return 0; PDEs shall be implemented is seperate files
////Real
////RBKernelSteadyState::computeQpResidual()
////{
////  return 0;
////}
//
////Real
////RBKernelSteadyState::computeQpJacobian()
////{
////  return 0;
////}
////#############################################################################################################################################
