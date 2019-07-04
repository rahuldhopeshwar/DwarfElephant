#ifndef DWARFELEPHANTSINGLECOORDINATEFILEREADER_H
#define DWARFELEPHANTSINGLECOORDINATEFILEREADER_H

#include "Function.h"

class DwarfElephantSingleCoordinateFileReader;

template <>
InputParameters validParams<DwarfElephantSingleCoordinateFileReader>();

///This function is responsible for reading boundary condition values from file.
class DwarfElephantSingleCoordinateFileReader : public Function
{
public:
  DwarfElephantSingleCoordinateFileReader(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;

protected:
  std::string _file;
  unsigned int _num_points;
  unsigned int _coordinate;
  std::vector<Real> _step_sizes;
  Real _tolerance;
  bool _interpolate;

  std::vector<Real> _coordinates;

private:
  void fileParser();
  Real findValue(Real _coord) const;
  Real interpolateValue(Real _coord) const;
};

#endif // DWARFELEPHANTBOUNDARYCONDITIONFILEREADER_H
