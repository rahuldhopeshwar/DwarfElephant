/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFELiftingFunctionKernelFunctionParameter.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFELiftingFunctionKernelFunctionParameter);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFELiftingFunctionKernelFunctionParameter>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("The class implements a lifting function.");
  params.addRequiredParam<FunctionName>("lifting_function", "Name of the lifting function two account for the inhomogeneous Dirichlet boundary conditions.");
  params.addRequiredParam<Real>("scale", "Defines the value of the scaling parameter.");
  params.addParam<Real>("norm_value", 1.0, "Defines the normalization value.");
  params.addRequiredParam<FunctionName>("function", "Name of the function that describes the parameter dependence.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFELiftingFunctionKernelFunctionParameter::DwarfElephantFELiftingFunctionKernelFunctionParameter(const InputParameters & parameters) :
  Diffusion(parameters),
  _lifting_function(&getFunction("lifting_function")),
  _scale(getParam<Real>("scale")),
  _norm_value(getParam<Real>("norm_value")),
  _func(getFunction("function"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation

Real
DwarfElephantFELiftingFunctionKernelFunctionParameter::computeQpResidual()
{
  Real dependency = _func.value(_fe_problem.time(),_q_point[_qp]);
  Real value = (_scale/_norm_value)* dependency*(_grad_test[_i][_qp]*(_lifting_function->gradient(_fe_problem.time(),_q_point[_qp])));
  return value;
}

Real
DwarfElephantFELiftingFunctionKernelFunctionParameter::computeQpJacobian()
{
   return 0.0;
}
