#include "DwarfElephantFECoupledEffectiveConductivity.h"

registerMooseObject("MooseApp", DwarfElephantFECoupledEffectiveConductivity);

template <>
InputParameters
validParams<DwarfElephantFECoupledEffectiveConductivity>()
{
  InputParameters params = validParams<GeneralPostprocessor>();

  params.addRequiredParam<PostprocessorName>("flux_top", "Name of the postprocessor that calculates the flux at the top surface");
  params.addRequiredParam<PostprocessorName>("flux_bottom", "Name of the postprocessor that calculates the flux at the bottom surface");
  params.addRequiredParam<PostprocessorName>("entropy", "Name of the postprocessor that calculates the entropy");
  params.addRequiredParam<PostprocessorName>("t_bar", "Name of the postprocessor that calculates the T_bar");
  params.addRequiredParam<Real>("thermal_conductivity", "Thermal conductivity of the unit.");

  return params;
}

DwarfElephantFECoupledEffectiveConductivity::DwarfElephantFECoupledEffectiveConductivity(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _flux_top(getParam<PostprocessorName>("flux_top")),
    _flux_bottom(getParam<PostprocessorName>("flux_bottom")),
    _entropy(getParam<PostprocessorName>("entropy")),
    _t_bar(getParam<PostprocessorName>("t_bar")),
    _thermal_conductivity(getParam<Real>("thermal_conductivity"))
{
  _flux_top_values = &getPostprocessorValueByName(_flux_top);
  _flux_bottom_values = &getPostprocessorValueByName(_flux_bottom);
  _entropy_values = &getPostprocessorValueByName(_entropy);
  _t_bar_values = &getPostprocessorValueByName(_t_bar);
}

void
DwarfElephantFECoupledEffectiveConductivity::initialize()
{
}

void
DwarfElephantFECoupledEffectiveConductivity::execute()
{
}

PostprocessorValue
DwarfElephantFECoupledEffectiveConductivity::getValue()
{
  Real effective_cond = 1.0;
  effective_cond *= abs(abs(*_flux_top_values) - abs(*_flux_bottom_values)) / *_t_bar_values;
  effective_cond *= sqrt(_thermal_conductivity / *_entropy_values);


  return effective_cond;
}
