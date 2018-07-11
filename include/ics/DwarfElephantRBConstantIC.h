#ifndef DWARFELEPHANTRBCONSTANTIC_H
#define DWARFELEPHANTRBCONSTANTIC_H

#include "DwarfElephantRBInitialCondition.h"

// System includes
#include <string>

// Forward Declarations
class DwarfElephantRBConstantIC;
class InputParameters;

namespace libMesh
{
class Point;
}

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<DwarfElephantRBConstantIC>();

/**
 * DwarfElephantRBConstantIC just returns a constant value.
 */
class DwarfElephantRBConstantIC : public DwarfElephantRBInitialCondition
{
public:
  DwarfElephantRBConstantIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:
  const Real _value;
};

#endif // DWARFELEPHANTRBCONSTANTIC_H
