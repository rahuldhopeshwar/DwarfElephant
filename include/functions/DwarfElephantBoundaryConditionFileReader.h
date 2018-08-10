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
  virtual RealGradient gradient(Real t, const Point & p) override;

protected:
  std::string _file;
  std::string _dx_file;
  std::string _dy_file;
  unsigned int _dimension;
  unsigned int _ID_data_layer;
  unsigned int _num_points;
  std::vector<unsigned int> _data_array_dimensions;
  std::vector<Real> _step_sizes;
  Real _tolerance;
  bool _interpolate;
  bool _access_multiple_times;
  bool _gradients;

  std::vector<Real> _x_coordinates;
  std::vector<Real> _y_coordinates;
  std::vector<Real> _z_coordinates;
  std::vector<Real> _variable_values;
  std::vector<Real> _dx;
  std::vector<Real> _dy;
  std::vector<std::vector<Real>> _data_array;

private:
  void fileParser();
  std::vector<Real>  fileParserGradients(std::string & file);
  Real findValue(Real _x_coord, Real _y_coord);
  Real interpolateValue(Real _x_coord, Real _y_coord);
};

#endif // DWARFELEPHANTBOUNDARYCONDITIONFILEREADER_H
