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
  int _dimension;
  int _num_points;

  std::vector<Real> _x_coordinates;
  std::vector<Real> _y_coordinates;
  std::vector<Real> _z_coordinates;
  std::vector<Real> _variable_values;
  // std::string _delimiter_to_replace;

private:
  void fileParser();
  Real findValue(Real _x_coord, Real _y_coord/*, Real _z_coord*/);
};

#endif // DWARFELEPHANTBOUNDARYCONDITIONFILEREADER_H
