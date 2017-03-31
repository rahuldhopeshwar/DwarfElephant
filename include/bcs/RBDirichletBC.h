#ifndef RBDIRICHLETBC_H
#define RBDIRICHLETBC_H

#include "RBNodalBC.h"

class RBDirichletBC;

template<>
InputParameters validParams<RBDirichletBC>();


class RBDirichletBC : public RBNodalBC
{
public:
  RBDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  /// The value for this BC
  const Real & _value;

  const DwarfElephantInitializeRBSystem & _initialize_rb_system;
};

#endif /* RBDIRICHLETBC_H */
