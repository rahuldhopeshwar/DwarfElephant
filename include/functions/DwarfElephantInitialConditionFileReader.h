#ifndef DWARFELEPHANTINITIALCONDITIONFILEREADER_H
#define DWARFELEPHANTINITIALCONDITIONFILEREADER_H

#include "Function.h"

class DwarfElephantInitialConditionFileReader;

template <>
InputParameters validParams<DwarfElephantInitialConditionFileReader>();

///This function is responsible for reading initial condition values from file.
class DwarfElephantInitialConditionFileReader : public Function
{
public:
  DwarfElephantInitialConditionFileReader(const InputParameters & parameters);
  virtual Real value(Real t, const Point & p) override;
  virtual Real value(const Node & n);

protected:
  const std::string _file;
  std::vector<std::vector<Real>> _values;

private:
  void fileReader();
};

#endif // DWARFELEPHANTINITIALCONDITIONFILEREADER_H
