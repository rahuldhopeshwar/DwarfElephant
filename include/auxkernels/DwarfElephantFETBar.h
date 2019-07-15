#ifndef DWARFELEPHANTFETBAR_H
#define DWARFELEPHANTFETBAR_H

#include "AuxKernel.h"
#include "Function.h"

// Forward Declarations
class DwarfElephantFETBar;

template<>
InputParameters validParams<DwarfElephantFETBar>();

class DwarfElephantFETBar : public AuxKernel
{
public:
  DwarfElephantFETBar(const InputParameters & parameters);

  virtual ~DwarfElephantFETBar() {}

protected:
  virtual Real computeValue() override;

  const Function & _T_top;
  const Function & _T_bottom;
};

#endif // DwarfElephantFETBar_H
