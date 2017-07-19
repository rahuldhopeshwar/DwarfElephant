#ifndef DWARFELEPHANTRBFUNCTIONNEUMANNBC_H
#define DWARFELEPHANTRBFUNCTIONNEUMANNBC_H

#include "DwarfElephantRBIntegratedBC.h"

//Forward Declarations
class DwarfElephantRBFunctionNeumannBC;
class Function;

template<>
InputParameters validParams<DwarfElephantRBFunctionNeumannBC>();

class DwarfElephantRBFunctionNeumannBC : public DwarfElephantRBIntegratedBC
{
public:
  DwarfElephantRBFunctionNeumannBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  
  Function & _func;
};

#endif // DWARFELEPHANTRBFUNCTIONNEUMANNBC_H
