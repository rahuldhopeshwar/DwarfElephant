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
  const Real & _value;
};

#endif /* RBDIRICHLETBC_H */
