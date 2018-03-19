/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFETHERMALCONDUCTIONNORMALIZED_H
#define DWARFELEPHANTFETHERMALCONDUCTIONNORMALIZED_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantFEThermalConduction.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEThermalConductionNormalized;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEThermalConductionNormalized>();

///-------------------------------------------------------------------------
class DwarfElephantFEThermalConductionNormalized : public DwarfElephantFEThermalConduction
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEThermalConductionNormalized(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
//  const MaterialProperty<Real> &_lambda; // for future use in case of varying material properties
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFETHERMALCONDUCTIONNORMALIZED_H
