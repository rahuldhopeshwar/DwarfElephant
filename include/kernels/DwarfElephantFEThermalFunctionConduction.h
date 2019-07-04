/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFETHERMALFUNCTIONCONDUCTION_H
#define DWARFELEPHANTFETHERMALFUNCTIONCONDUCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEThermalFunctionConduction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEThermalFunctionConduction>();

///This Kernel is implements a thermal conduction problem using the full Finite Element solution. It is included in this package for validation purposes.
class DwarfElephantFEThermalFunctionConduction : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEThermalFunctionConduction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
//  const MaterialProperty<Real> &_lambda; // for future use in case of varying material properties
  Real _lambda;
  Real _norm_value;

  const Function & _func;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFETHERMALFUNCTIONCONDUCTION_H
