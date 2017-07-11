#ifndef DWARFELEPAHNTRBNEUMANNBC_H
#define DWARFELEPHANTRBNEUMANNBC_H

#include "DwarfElephantRBIntegratedBC.h"


class DwarfElephantRBNeumannBC;

template<>
InputParameters validParams<DwarfElephantRBNeumannBC>();

class DwarfElephantRBNeumannBC : public DwarfElephantRBIntegratedBC
{
public:

  DwarfElephantRBNeumannBC(const InputParameters & parameters);


protected:
  virtual Real computeQpResidual() override;
  const Real & _value;
};


#endif //DWARFELEPHANTRBNEUMANNBC_H
