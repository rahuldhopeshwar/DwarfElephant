/**
 * This Kernel is implements a constant radiogenic heat production using the
 * Reduced Basis solution. It is required if you do want the radiogenic heat
 * production stay fixed for all parameters.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBCONSTANTRADIOGENICHEATPRODUCTION_H
#define DWARFELEPHANTRBCONSTANTRADIOGENICHEATPRODUCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantRBKernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBConstantRadiogenicHeatProduction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBConstantRadiogenicHeatProduction>();

///This Kernel is implements a constant radiogenic heat production using the  Reduced Basis solution. It is required if you do want the radiogenic heat production stay fixed for all parameters.
class DwarfElephantRBConstantRadiogenicHeatProduction : public DwarfElephantRBKernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBConstantRadiogenicHeatProduction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  Real _radiogenic_heat_production;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCONSTANTRADIOGENICHEATPRODUCTION_H
