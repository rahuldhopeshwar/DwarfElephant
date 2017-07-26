#ifndef DWARFELEPHANTRBPRESETNODALBC_H
#define DWARFELEPHANTRBPRESETNODALBC_H

#include "DwarfElephantRBDirichletBC.h"
#include "MooseVariable.h"

class DwarfElephantRBPresetNodalBC;

template<>
InputParameters validParams<DwarfElephantRBPresetNodalBC>();


class DwarfElephantRBPresetNodalBC : public DwarfElephantRBNodalBC
{
public:
  DwarfElephantRBPresetNodalBC(const InputParameters & parameters);

  void computeValue(NumericVector<Number> & current_solution);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpValue();

};

#endif /* DWARFELEPHANTRBPRESETBC_H */
