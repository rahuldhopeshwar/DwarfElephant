/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTCONDUCTION_H
#define DWARFELEPHANTCONDUCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantConduction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantConduction>();

///-------------------------------------------------------------------------
class DwarfElephantConduction : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantConduction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  const MaterialProperty<Real> &_lambda;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTCONDUCTION_H
