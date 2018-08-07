/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTLIFTINGFUNCTION_H
#define DWARFELEPHANTLIFTINGFUNCTION_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantLiftingFunction;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantLiftingFunction>();

///-------------------------------------------------------------------------
class DwarfElephantLiftingFunction : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantLiftingFunction(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  Function * _lifting_function;
  Real _scale;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTLIFTINGFUNCTION_H
