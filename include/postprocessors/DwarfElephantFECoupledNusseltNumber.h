#include "GeneralPostprocessor.h"

class DwarfElephantFECoupledNusseltNumber;

template <>
InputParameters validParams<DwarfElephantFECoupledNusseltNumber>();

class DwarfElephantFECoupledNusseltNumber : public GeneralPostprocessor
{
public:
  DwarfElephantFECoupledNusseltNumber(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual PostprocessorValue getValue() override;

protected:
  const std::vector<PostprocessorName> & _pp_names;
  const unsigned int _n_pp;
  const Real _thermal_conductivity;

  std::vector<const PostprocessorValue *> _pp_values;
};
