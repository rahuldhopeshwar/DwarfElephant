#ifndef DWARFELEPHANTLIFTINGFUNCTION_H
#define DWARFELEPHANTLIFTINGFUNCTION_H

#include "Function.h"
#include "FunctionInterface.h"

class DwarfElephantLiftingFunction;

template <>
InputParameters validParams<DwarfElephantLiftingFunction>();

class DwarfElephantLiftingFunction : public Function, protected FunctionInterface
{
public:
  DwarfElephantLiftingFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) override;
  virtual RealGradient gradient(Real t, const Point & p) override;
  virtual void initialSetup() override;

protected:
  Function * _boundary_function;
  unsigned int _ID_reference_layer;
  unsigned int _ID_data_layer;
  std::vector<unsigned int> _data_array_dimensions;
  std::vector<Real> _step_sizes;
  Real _tolerance;
  std::vector<Real> _x_coord_reference_layer;
  std::vector<Real> _y_coord_reference_layer;
  std::vector<Real> _z_coord_reference_layer;
  std::vector<std::vector<Real>> _data_array;
  bool _access_multiple_times;

private:
  Real interpolateZRef(Real _x_coord, Real _y_coord);
};

#endif // DWARFELEPHANTLIFTINGFUNCTION_H
