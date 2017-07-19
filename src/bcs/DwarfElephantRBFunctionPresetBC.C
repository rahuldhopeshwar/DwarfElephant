#include "DwarfElephantRBFunctionPresetBC.h"
#include "Function.h"

template<>
InputParameters validParams<DwarfElephantRBFunctionPresetBC>()
{
  InputParameters params = validParams<DwarfElephantRBPresetNodalBC>();
  params.addRequiredParam<FunctionName>("function", "The forcing function.");
  return params;
}

DwarfElephantRBFunctionPresetBC::DwarfElephantRBFunctionPresetBC(const InputParameters & parameters) :
    DwarfElephantRBPresetNodalBC(parameters),
    _func(getFunction("function"))
{
}

Real
DwarfElephantRBFunctionPresetBC::computeQpValue()
{
  return _func.value(_t, *_current_node);
}

