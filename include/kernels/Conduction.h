/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef CONDUCTION_H
#define CONDUCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Kernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class Conduction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<Conduction>();

///-------------------------------------------------------------------------
class Conduction : public Kernel
{
//----------------------------------PUBLIC----------------------------------
public:
  Conduction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  const MaterialProperty<Real> &_lambda;
};

///-------------------------------------------------------------------------
#endif // CONDUCTION_H
