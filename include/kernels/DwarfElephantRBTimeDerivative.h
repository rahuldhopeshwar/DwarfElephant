#ifndef DWARFELEPHANTRBTIMEDERIVATIVE_H
#define DWARFELEPHANTRBTIMEDERIVATIVE_H

#include "DwarfElephantRBTimeKernel.h"

// Forward Declaration
class DwarfElephantRBTimeDerivative;

template <>
InputParameters validParams<DwarfElephantRBTimeDerivative>();

class DwarfElephantRBTimeDerivative : public DwarfElephantRBTimeKernel
{
public:
  DwarfElephantRBTimeDerivative(const InputParameters & parameters);

  virtual void computeJacobian() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  bool _lumping;
};

#endif // DWARFELEPHANTRBTIMEDERIVATIVE_H
