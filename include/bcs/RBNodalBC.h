#ifndef RBNODALBC_H
#define RBNODALBC_H

#include "NodalBC.h"

// Forward declarations
class RBNodalBC;

template<>
InputParameters validParams<RBNodalBC>();

class RBNodalBC :
  public NodalBC
{
public:
  RBNodalBC(const InputParameters & parameters);

  virtual void computeJacobian();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};

#endif /* RBNODALBC_H */
