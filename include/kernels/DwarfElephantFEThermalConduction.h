/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFETHERMALCONDUCTION_H
#define DWARFELEPHANTFETHERMALCONDUCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEThermalConduction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEThermalConduction>();

///This Kernel is implements a thermal conduction problem using the full Finite Element solution. It is included in this package for validation purposes.
class DwarfElephantFEThermalConduction : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEThermalConduction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
//  const MaterialProperty<Real> &_lambda; // for future use in case of varying material properties
  Real _lambda;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFETHERMALCONDUCTION_H
