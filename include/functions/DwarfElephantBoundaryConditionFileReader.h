#ifndef DWARFELEPHANTBOUNDARYCONDITIONFILEREADER_H
#define DWARFELEPHANTBOUNDARYCONDITIONFILEREADER_H

#include "Function.h"

class DwarfElephantBoundaryConditionFileReader;

template <>
InputParameters validParams<DwarfElephantBoundaryConditionFileReader>();

class DwarfElephantBoundaryConditionFileReader : public Function
{
public:
  DwarfElephantBoundaryConditionFileReader(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) override;

protected:
  const std::string _file;

private:
  void fileReader();
};

#endif // DWARFELEPHANTBOUNDARYCONDITIONFILEREADER_H
