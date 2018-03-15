#ifndef DWARFELEPHANTFILEIC_H
#define DWARFELEPHANTFILEIC_H

#include "FunctionIC.h"

#include <string>

// Forward Declarations
class DwarfElephantFileIC;

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<DwarfElephantFileIC>();

class DwarfElephantFileIC : public FunctionIC
{
public:
  DwarfElephantFileIC(const InputParameters & parameters);

protected:
  Real f();
  virtual Real value(const Point & p) override;
};

#endif // DWARFELEPHANTFILEIC_H
