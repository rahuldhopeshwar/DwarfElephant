#ifndef DWARFELEPHANTRBFUNCTIONPENALTYDIRICHLETBC_H
#define DWARFELEPHANTRBFUNCTIONPENALTYDIRICHLETBC_H

#include "DwarfElephantRBIntegratedBC.h"

class DwarfElephantRBFunctionPenaltyDirichletBC;
class Function;

template <>
InputParameters validParams<DwarfElephantRBFunctionPenaltyDirichletBC>();

/**
 * A different approach to applying Dirichlet BCs
 *
 * uses \f$ \int(p u \cdot \phi)=\int(p f \cdot \phi)\f$ on \f$d\omega\f$
 *
 */

class DwarfElephantRBFunctionPenaltyDirichletBC : public DwarfElephantRBIntegratedBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  DwarfElephantRBFunctionPenaltyDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  Function & _func;

private:
  Real _p;
};

#endif
