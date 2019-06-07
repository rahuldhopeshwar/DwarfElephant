#ifndef DWARFELEPHANTFEPENALTYFACTORFUNCTIONDIRICHLETBC_H
#define DWARFELEPHANTFEPENALTYFACTORFUNCTIONDIRICHLETBC_H

#include "IntegratedBC.h"

class DwarfElephantFEPenaltyFactorFunctionDirichletBC;
class Function;

template <>
InputParameters validParams<DwarfElephantFEPenaltyFactorFunctionDirichletBC>();


class DwarfElephantFEPenaltyFactorFunctionDirichletBC : public IntegratedBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  DwarfElephantFEPenaltyFactorFunctionDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const Function & _func;

private:
  Real _value;
};

#endif
