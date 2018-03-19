/**
 * This Kernel is implements a darcy flow problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFEDARCYNORMALIZED_H
#define DWARFELEPHANTFEDARCYNORMALIZED_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantFEDarcy.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEDarcyNormalized;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEDarcyNormalized>();

///-------------------------------------------------------------------------
class DwarfElephantFEDarcyNormalized : public DwarfElephantFEDarcy
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEDarcyNormalized(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFEDARCYNORMALIZED_H
