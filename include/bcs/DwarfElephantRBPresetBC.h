#ifndef DWARFELEPHANTRBPRESETBC_H
#define DWARFELEPHANTRBPRESETBC_H

#include "DwarfElephantRBPresetNodalBC.h"

class DwarfElephantRBPresetBC;

template<>
InputParameters validParams<DwarfElephantRBPresetBC>();


class DwarfElephantRBPresetBC : public DwarfElephantRBPresetNodalBC
{
public:
  DwarfElephantRBPresetBC(const InputParameters & parameters);

protected:
  virtual Real computeQpValue() override;

  const Real & _value;
};

#endif /* DWARFELEPHANTRBPRESETBC_H */
