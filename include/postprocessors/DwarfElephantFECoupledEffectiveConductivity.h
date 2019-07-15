#include "GeneralPostprocessor.h"

class DwarfElephantFECoupledEffectiveConductivity;

template <>
InputParameters validParams<DwarfElephantFECoupledEffectiveConductivity>();

class DwarfElephantFECoupledEffectiveConductivity : public GeneralPostprocessor
{
public:
  DwarfElephantFECoupledEffectiveConductivity(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual PostprocessorValue getValue() override;

protected:
  const PostprocessorName & _flux_top;
  const PostprocessorName & _flux_bottom;
  const PostprocessorName & _entropy;
  const PostprocessorName & _t_bar;
  const Real _thermal_conductivity;


  const PostprocessorValue * _flux_top_values;
  const PostprocessorValue * _flux_bottom_values;
  const PostprocessorValue * _entropy_values;
  const PostprocessorValue * _t_bar_values;
};
