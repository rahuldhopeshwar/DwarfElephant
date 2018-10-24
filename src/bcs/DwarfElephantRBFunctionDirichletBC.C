#include "DwarfElephantRBFunctionDirichletBC.h"
#include "Function.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBFunctionDirichletBC);

template<>
InputParameters validParams<DwarfElephantRBFunctionDirichletBC>()
{
  InputParameters params = validParams<DwarfElephantRBNodalBC>();
  params.addRequiredParam<FunctionName>("function", "The forcing function.");
  return params;
}

DwarfElephantRBFunctionDirichletBC::DwarfElephantRBFunctionDirichletBC(const InputParameters & parameters) :
    DwarfElephantRBNodalBC(parameters),
    _func(getFunction("function"))
{
}

Real
DwarfElephantRBFunctionDirichletBC::f()
{
  return _func.value(_t, *_current_node);
}

Real
DwarfElephantRBFunctionDirichletBC::computeQpResidual()
{
  return _u[_qp]-f();
}
