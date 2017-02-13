#ifndef KERNELOUTPUTAUX_H
#define KERNELOUTPUTAUX_H

#include "AuxKernel.h"


//Forward Declarations
class KernelOutputAux;

template<>
InputParameters validParams<KernelOutputAux>();

class KernelOutputAux : public AuxKernel
{
public:

  KernelOutputAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

//  const VariableValue & _coupled_val;
//
//  Real _value;
};

#endif //KERNELOUTPUTAUX_H
