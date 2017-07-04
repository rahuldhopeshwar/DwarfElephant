/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTCONDUCTIONLIFTINGFUNCTION_H
#define DWARFELEPHANTCONDUCTIONLIFTINGFUNCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantConductionLiftingFunction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantConductionLiftingFunction>();

///-------------------------------------------------------------------------
class DwarfElephantConductionLiftingFunction : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantConductionLiftingFunction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
 

  /* Attributes */
  const MaterialProperty<Real> &_lambda;
  
  Function * _lifting_function;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTCONDUCTIONLIFTINGFUNCTION_H
