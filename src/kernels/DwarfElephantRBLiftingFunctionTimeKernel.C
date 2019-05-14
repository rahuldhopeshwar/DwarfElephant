/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBLiftingFunctionTimeKernel.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBLiftingFunctionTimeKernel);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBLiftingFunctionTimeKernel>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  params.addClassDescription("The class implements a lifting function.");
  params.addRequiredParam<FunctionName>("lifting_function", "Name of the lifting function two account for the inhomogeneous Dirichlet boundary conditions.");

  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBLiftingFunctionTimeKernel::DwarfElephantRBLiftingFunctionTimeKernel(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _lifting_function(&getFunction("lifting_function"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation

Real
DwarfElephantRBLiftingFunctionTimeKernel::computeQpResidual()
{
  return  -_lifting_function->timeDerivative(_fe_problem.time(),_q_point[_qp])* _test[_i][_qp];
}

Real
DwarfElephantRBLiftingFunctionTimeKernel::computeQpJacobian()
{
   return 0.0;
}
