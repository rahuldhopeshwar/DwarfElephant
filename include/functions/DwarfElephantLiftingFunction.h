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
  virtual void initialSetup() override;
  //TODO: add the function for the gradient

protected:
  Function & _boundary_function;
  const std::vector<dof_id_type> * _bottom_nodes;
};

#endif // DWARFELEPHANTLIFTINGFUNCTION_H
