/**
 * This Kernel is implements a constant radiogenic heat production using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFECONSTANTRADIOGENICHEATPRODUCTION_H
#define DWARFELEPHANTFECONSTANTRADIOGENICHEATPRODUCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Kernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEConstantRadiogenicHeatProduction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEConstantRadiogenicHeatProduction>();

///This Kernel is implements a constant radiogenic heat production using the full Finite Element solution. It is included in this package for validation purposes.
class DwarfElephantFEConstantRadiogenicHeatProduction : public Kernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEConstantRadiogenicHeatProduction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  Real _radiogenic_heat_production;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFECONSTANTRADIOGENICHEATPRODUCTION_H
