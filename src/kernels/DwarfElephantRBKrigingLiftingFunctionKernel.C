/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBKrigingLiftingFunctionKernel.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBKrigingLiftingFunctionKernel);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBKrigingLiftingFunctionKernel>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  params.addClassDescription("The class implements a lifting function.");
  params.addRequiredParam<FunctionName>("lifting_function_1", "Name of the first lifting function two account for the inhomogeneous Dirichlet boundary conditions.");
  params.addRequiredParam<FunctionName>("lifting_function_2", "Name of the second lifting function two account for the inhomogeneous Dirichlet boundary conditions.");
  params.addRequiredParam<Real>("range", "Range of the semivariogramm.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBKrigingLiftingFunctionKernel::DwarfElephantRBKrigingLiftingFunctionKernel(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _range(getParam<Real>("range")),
  _lifting_function_1(&getFunction("lifting_function_1")),
  _lifting_function_2(&getFunction("lifting_function_2"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation

Real
DwarfElephantRBKrigingLiftingFunctionKernel::computeQpResidual()
{
  Real h = sqrt(pow(_q_point[_qp](0),2)+
           pow(_q_point[_qp](1),2)+
           pow(_q_point[_qp](2),2));
  Real value;

  if(h>0 && h<_range)
    value = -(_grad_test[_i][_qp]*(_lifting_function_1->gradient(_fe_problem.time(),_q_point[_qp])));
  else
    value = -(_grad_test[_i][_qp]*(_lifting_function_2->gradient(_fe_problem.time(),_q_point[_qp])));
  return value;
}

Real
DwarfElephantRBKrigingLiftingFunctionKernel::computeQpJacobian()
{
   return 0.0;
}
