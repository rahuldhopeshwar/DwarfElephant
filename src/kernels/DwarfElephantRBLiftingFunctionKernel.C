/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBLiftingFunctionKernel.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBLiftingFunctionKernel);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBLiftingFunctionKernel>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  params.addClassDescription("The class implements a lifting function.");
  params.addRequiredParam<FunctionName>("lifting_function", "Name of the lifting function two account for the inhomogeneous Dirichlet boundary conditions.");

  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBLiftingFunctionKernel::DwarfElephantRBLiftingFunctionKernel(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _lifting_function(&getFunction("lifting_function"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation

Real
DwarfElephantRBLiftingFunctionKernel::computeQpResidual()
{
  return  -(_grad_test[_i][_qp]*(_lifting_function->gradient(_fe_problem.time(),_q_point[_qp])));
}

Real
DwarfElephantRBLiftingFunctionKernel::computeQpJacobian()
{
   return 0.0;
}
