/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEKrigingLiftingFunctionKernel.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFEKrigingLiftingFunctionKernel);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEKrigingLiftingFunctionKernel>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("The class implements a lifting function.");
  params.addRequiredParam<FunctionName>("lifting_function_1", "Name of the first lifting function two account for the inhomogeneous Dirichlet boundary conditions.");
  params.addRequiredParam<FunctionName>("lifting_function_2", "Name of the second lifting function two account for the inhomogeneous Dirichlet boundary conditions.");
  params.addRequiredParam<Real>("scale", "Defines the value of the scaling parameter.");
  params.addRequiredParam<Real>("range", "Range of the semivariogramm.");
  params.addParam<Real>("norm_value", 1.0, "Defines the normalization value.");

  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEKrigingLiftingFunctionKernel::DwarfElephantFEKrigingLiftingFunctionKernel(const InputParameters & parameters) :
  Diffusion(parameters),
  _range(getParam<Real>("range")),
  _lifting_function_1(&getFunction("lifting_function_1")),
  _lifting_function_2(&getFunction("lifting_function_2")),
  _scale(getParam<Real>("scale")),
  _norm_value(getParam<Real>("norm_value"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation

Real
DwarfElephantFEKrigingLiftingFunctionKernel::computeQpResidual()
{
  Real h = sqrt(pow(_q_point[_qp](0),2)+
           pow(_q_point[_qp](1),2)+
           pow(_q_point[_qp](2),2));
  Real value;

  if(h>0 && h<_range)
    value = (_scale/_norm_value)*(_grad_test[_i][_qp]*(_lifting_function_1->gradient(_fe_problem.time(),_q_point[_qp])));
  else
    value = (_scale/_norm_value)*(_grad_test[_i][_qp]*(_lifting_function_2->gradient(_fe_problem.time(),_q_point[_qp])));
  return value;
}

Real
DwarfElephantFEKrigingLiftingFunctionKernel::computeQpJacobian()
{
   return 0.0;
}
