#ifndef DWARFELEPAHNTRBNEUMANNBCND_H
#define DWARFELEPHANTRBNEUMANNBCND_H

#include "DwarfElephantRBIntegratedBC.h"


class DwarfElephantRBNeumannBCND;

template<>
InputParameters validParams<DwarfElephantRBNeumannBCND>();

class DwarfElephantRBNeumannBCND : public DwarfElephantRBIntegratedBC
{
public:

  DwarfElephantRBNeumannBCND(const InputParameters & parameters);


protected:
  virtual Real computeQpResidual() override;
  const Real & _value;
  Real _u_ref;
  Real _l_ref;
};


#endif //DWARFELEPHANTRBNEUMANNBC_H
