/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTLIFTINGFUNCTIONKERNEL_H
#define DWARFELEPHANTLIFTINGFUNCTIONKERNEL_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantLiftingFunctionKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantLiftingFunctionKernel>();

///-------------------------------------------------------------------------
class DwarfElephantLiftingFunctionKernel : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantLiftingFunctionKernel(const InputParameters & parameters);

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
#endif // DWARFELEPHANTLIFTINGFUNCTIONKERNEL_H
