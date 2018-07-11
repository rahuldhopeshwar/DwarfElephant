#include "Function.h"

#include "DwarfElephantRBFunctionIC.h"

template <>
InputParameters
validParams<DwarfElephantRBFunctionIC>()
{
  InputParameters params = validParams<DwarfElephantRBInitialCondition>();
  params.addRequiredParam<FunctionName>("function", "The initial condition function.");

  params.addClassDescription("An initial condition that uses a normal function of x, y, z to "
                             "produce values (and optionally gradients) for a field variable."
                             "This class is only required if you are using parameter-dependent "
                             "initial conditions.");
  return params;
}

DwarfElephantRBFunctionIC::DwarfElephantRBFunctionIC(const InputParameters & parameters)
  : DwarfElephantRBInitialCondition(parameters),
   _func(getFunction("function"))
{
}

Real
DwarfElephantRBFunctionIC::value(const Point & p)
{
  return _func.value(_t, p);
}

RealGradient
DwarfElephantRBFunctionIC::gradient(const Point & p)
{
  return _func.gradient(_t, p);
}
