/**
 * This Kernel is implements a darcy flow problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFEDARCY_H
#define DWARFELEPHANTFEDARCY_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Kernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEDarcy;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEDarcy>();

///-------------------------------------------------------------------------
class DwarfElephantFEDarcy : public Kernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEDarcy(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  Real _permeability;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFEDARCY_H
