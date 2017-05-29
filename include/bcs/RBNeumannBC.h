#ifndef RBNEUMANNBC_H
#define RBNEUMANNBC_H

#include "RBIntegratedBC.h"


class RBNeumannBC;

template<>
InputParameters validParams<RBNeumannBC>();

class RBNeumannBC : public RBIntegratedBC
{
public:

  RBNeumannBC(const InputParameters & parameters);


protected:
  virtual Real computeQpResidual() override;
  const Real & _value;
};


#endif //RBNEUMANNBC_H
