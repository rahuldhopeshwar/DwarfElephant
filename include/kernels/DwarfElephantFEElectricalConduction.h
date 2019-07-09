/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFEELECTRICALCONDUCTION_H
#define DWARFELEPHANTFEELECTRICALCONDUCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEElectricalConduction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEElectricalConduction>();

///This Kernel is implements a thermal conduction problem using the full Finite Element solution. It is included in this package for validation purposes.
class DwarfElephantFEElectricalConduction : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEElectricalConduction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  Real _resistivity;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFEELECTRICALCONDUCTION_H
