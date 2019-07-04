/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#ifndef DWARFELEPHANTFECOUPLEDTEMPERATUREDIFFUSION_H
#define DWARFELEPHANTFECOUPLEDTEMPERATUREDIFFUSION_H

#include "Diffusion.h"

//Forward Declarations
class DwarfElephantFECoupledTemperatureDiffusion;

/* This class extends the Diffusion kernel to multiply by a coefficient
 * read from the input file
 */
template<>
InputParameters validParams<DwarfElephantFECoupledTemperatureDiffusion>();

class DwarfElephantFECoupledTemperatureDiffusion : public Diffusion
{
public:

  DwarfElephantFECoupledTemperatureDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  Real _diffusivity;
};
#endif //DwarfElephantFECoupledTemperatureDiffusion_H
