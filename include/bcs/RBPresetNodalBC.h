#ifndef RBPRESETNODALBC_H
#define RBPRESETNODALBC_H

#include "RBNodalBC.h"

class RBPresetNodalBC;

template<>
InputParameters validParams<RBPresetNodalBC>();


class RBPresetNodalBC : public RBNodalBC
{
public:
  RBPresetNodalBC(const InputParameters & parameters);

  void computeValue(NumericVector<Number> & current_solution);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpValue();

};

#endif /* RBPRESETBC_H */
