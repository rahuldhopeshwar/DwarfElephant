/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#ifndef DWARFELEPHANTFECOUPLEDPRESSUREDIFFUSION_H
#define DWARFELEPHANTFECOUPLEDPRESSUREDIFFUSION_H

#include "Diffusion.h"

//Forward Declarations
class DwarfElephantFECoupledPressureDiffusion;

/* This class extends the Diffusion kernel to multiply by a coefficient
 * read from the input file
 */
template<>
InputParameters validParams<DwarfElephantFECoupledPressureDiffusion>();

class DwarfElephantFECoupledPressureDiffusion : public Diffusion
{
public:

  DwarfElephantFECoupledPressureDiffusion(const InputParameters & parameters);
  virtual ~DwarfElephantFECoupledPressureDiffusion() {}

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  const VariableValue & _temp;
  unsigned _temp_var_num;
  const VariableGradient & _grad_temp;
  const MaterialProperty<Real> & _Ra;
  unsigned _component;
};
#endif //DwarfElephantFECoupledPressureDiffusion_H
