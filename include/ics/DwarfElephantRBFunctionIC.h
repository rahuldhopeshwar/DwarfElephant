#ifndef DWARFELEPHANTRBFUNCTIONIC_H
#define DWARFELEPHANTRBFUNCTIONIC_H

#include "DwarfElephantRBInitialCondition.h"

#include <string>

// Forward Declarations
class DwarfElephantRBFunctionIC;
class Function;
class InputParameters;

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<DwarfElephantRBFunctionIC>();

/**
 * Defines a boundary condition that forces the value to be a user specified
 * function at the boundary.
 */
class DwarfElephantRBFunctionIC : public DwarfElephantRBInitialCondition
{
public:
  DwarfElephantRBFunctionIC(const InputParameters & parameters);

protected:
  /**
   * Evaluate the function at the current quadrature point and time step.
   */
  Real f();

  /**
   * The value of the variable at a point.
   */
  virtual Real value(const Point & p) override;

  /**
   * The value of the gradient at a point.
   */
  virtual RealGradient gradient(const Point & p) override;

  Function & _func;
};

#endif // DWARFELEPHANTRBFUNCTIONIC_H
