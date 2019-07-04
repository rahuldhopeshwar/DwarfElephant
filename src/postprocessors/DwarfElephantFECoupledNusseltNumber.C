#include "DwarfElephantFECoupledNusseltNumber.h"

registerMooseObject("MooseApp", DwarfElephantFECoupledNusseltNumber);

template <>
InputParameters
validParams<DwarfElephantFECoupledNusseltNumber>()
{
  InputParameters params = validParams<GeneralPostprocessor>();

  params.addRequiredParam<std::vector<PostprocessorName>>("pp_names", "List of post-processors");
  params.addRequiredParam<Real>("thermal_conductivity", "Thermal conductivity of the unit.");

  return params;
}

DwarfElephantFECoupledNusseltNumber::DwarfElephantFECoupledNusseltNumber(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _pp_names(getParam<std::vector<PostprocessorName>>("pp_names")),
    _n_pp(_pp_names.size()),
    _thermal_conductivity(getParam<Real>("thermal_conductivity"))
{
  _pp_values.resize(_n_pp);
  for (unsigned int i = 0; i < _n_pp; i++)
    _pp_values[i] = &getPostprocessorValueByName(_pp_names[i]);
}

void
DwarfElephantFECoupledNusseltNumber::initialize()
{
}

void
DwarfElephantFECoupledNusseltNumber::execute()
{
}

PostprocessorValue
DwarfElephantFECoupledNusseltNumber::getValue()
{
  Real nusselt_number = 1./_thermal_conductivity;
  for (unsigned int i = 0; i < _n_pp; i++)
    nusselt_number *= *(_pp_values[i]);

  return nusselt_number;
}
