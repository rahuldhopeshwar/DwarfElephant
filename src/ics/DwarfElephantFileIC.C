#include "DwarfElephantFileIC.h"
#include "DwarfElephantInitialConditionFileReader.h"
#include "Function.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFileIC);

template <>
InputParameters
validParams<DwarfElephantFileIC>()
{
  InputParameters params = validParams<FunctionIC>();
  return params;
}

DwarfElephantFileIC::DwarfElephantFileIC(const InputParameters & parameters)
  : FunctionIC(parameters)
{
}

Real
DwarfElephantFileIC::value(const Point & /*p*/)
{
  DwarfElephantInitialConditionFileReader & _file_func = cast_ref<DwarfElephantInitialConditionFileReader &>(_func);
  return _file_func.value(*_current_node);
}
