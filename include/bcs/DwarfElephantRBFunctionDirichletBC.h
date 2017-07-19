#ifndef DWARFELEPHANTRBFUNCTIONDIRICHLETBC_H
#define DWARFELEPHANTRBFUNCTIONDIRICHLETBC_H

#include "DwarfElephantRBNodalBC.h"

//Forward Declarations
class DwarfElephantRBFunctionDirichletBC;
class Function;

template<>
InputParameters validParams<DwarfElephantRBFunctionDirichletBC>();


class DwarfElephantRBFunctionDirichletBC : public DwarfElephantRBNodalBC
{
public:
  DwarfElephantRBFunctionDirichletBC(const InputParameters & parameters);

protected:
  Real f();

  virtual Real computeQpResidual() override;

  Function & _func;
};

#endif //DWARFELEPHANTRBFUNCTIONDIRICHLETBC_H
