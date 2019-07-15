#include "DwarfElephantFETBar.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFETBar);

template<>
InputParameters validParams<DwarfElephantFETBar>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredParam<FunctionName>("T_top", "A function that describes the temperature at the top of the model.");
  params.addRequiredParam<FunctionName>("T_bottom", "A function that describes the temperature at the bottom of the model.");

  return params;
}

DwarfElephantFETBar::DwarfElephantFETBar(const InputParameters & parameters) :
    AuxKernel(parameters),
    _T_top(getFunction("T_top")),
    _T_bottom(getFunction("T_bottom"))
{
}

Real
DwarfElephantFETBar::computeValue()
{
  Real _T_bar = (_T_top.value(_c_fe_problem.time(),_q_point[_qp])+
                 _T_bottom.value(_c_fe_problem.time(),_q_point[_qp]))/2;

  return _T_bar;

}
