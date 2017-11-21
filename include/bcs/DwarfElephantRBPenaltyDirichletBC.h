#ifndef DWARFELEPHANTRBPENALTYDIRICHLETBC_H
#define DWARFELEPHANTRBPENALTYDIRICHLETBC_H

#include "DwarfElephantRBIntegratedBC.h"

class DwarfElephantRBPenaltyDirichletBC;
class Function;

template <>
InputParameters validParams<DwarfElephantRBPenaltyDirichletBC>();

class DwarfElephantRBPenaltyDirichletBC : public DwarfElephantRBIntegratedBC
{
public:
  DwarfElephantRBPenaltyDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  Real _p;
  const Real & _v;
};

#endif // DWARFELEPHANTRBPENEALTYDIRICHLETBC_H
