#ifndef DWARFELEPHANTRBDIRICHLETBC_H
#define DWARFELEPHANTRBDIRICHLETBC_H

#include "DwarfElephantRBNodalBC.h"

class DwarfElephantRBDirichletBC;

template<>
InputParameters validParams<DwarfElephantRBDirichletBC>();


class DwarfElephantRBDirichletBC : public DwarfElephantRBNodalBC
{
public:
  DwarfElephantRBDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  const Real & _value;
};

#endif /* DWARFELEPHANTRBDIRICHLETBC_H */
