#ifndef RBPRESETBC_H
#define RBPRESETBC_H

#include "RBPresetNodalBC.h"

class RBPresetBC;

template<>
InputParameters validParams<RBPresetBC>();


class RBPresetBC : public RBPresetNodalBC
{
public:
  RBPresetBC(const InputParameters & parameters);

protected:
  virtual Real computeQpValue() override;

  const Real & _value;
};

#endif /* RBPRESETBC_H */
