#include "DwarfElephantFEKrigingFunctionDirichletBC.h"
#include "Function.h"

registerMooseObject("MooseApp", DwarfElephantFEKrigingFunctionDirichletBC);

template <>
InputParameters
validParams<DwarfElephantFEKrigingFunctionDirichletBC>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredParam<FunctionName>("function_1", "Option 1 for the forcing function.");
  params.addRequiredParam<FunctionName>("function_2", "Option 2 for the forcing function.");
  params.addRequiredParam<Real>("range", "Range of the semivariogramm.");
  return params;
}

DwarfElephantFEKrigingFunctionDirichletBC::DwarfElephantFEKrigingFunctionDirichletBC(const InputParameters & parameters)
  : NodalBC(parameters),
   _range(getParam<Real>("range")),
   _func_1(getFunction("function_1")),
   _func_2(getFunction("function_2"))

{
}

Real
DwarfElephantFEKrigingFunctionDirichletBC::f()
{
  Real h = sqrt(pow(_current_node->operator()(0),2)+
           pow(_current_node->operator()(1),2)+
           pow(_current_node->operator()(2),2));
  Real value;

  if(h>0 && h<_range)
    value = _func_1.value(_t, *_current_node);
  else
    value = _func_2.value(_t, *_current_node);

  return value;
}

Real
DwarfElephantFEKrigingFunctionDirichletBC::computeQpResidual()
{
  return _u[_qp] - f();
}
