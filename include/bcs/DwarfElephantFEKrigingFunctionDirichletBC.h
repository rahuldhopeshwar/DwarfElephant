#ifndef DWARFELEPHANTFEKRIGINGFUNCTIONDIRICHLETBC_H
#define DWARFELEPHANTFEKRIGINGFUNCTIONDIRICHLETBC_H

#include "NodalBC.h"

// Forward Declarations
class DwarfElephantFEKrigingFunctionDirichletBC;
class Function;

template <>
InputParameters validParams<DwarfElephantFEKrigingFunctionDirichletBC>();

class DwarfElephantFEKrigingFunctionDirichletBC : public NodalBC
{
public:
  DwarfElephantFEKrigingFunctionDirichletBC(const InputParameters & parameters);

protected:
  Real f();

  virtual Real computeQpResidual() override;

  /// The function being used for evaluation
  Real _range;
  Function & _func_1;
  Function & _func_2;
};

#endif // DWARFELEPHANTFEKRIGINGFUNCTIONDIRICHLETBC_H
