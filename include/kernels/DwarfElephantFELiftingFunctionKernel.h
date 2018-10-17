/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFELIFTINGFUNCTIONKERNEL_H
#define DWARFELEPHANTFELIFTINGFUNCTIONKERNEL_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFELiftingFunctionKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFELiftingFunctionKernel>();

///This Kernel is implements the concept of the lifting function to avoid problems caused by non-homogenous DirichletBC.
class DwarfElephantFELiftingFunctionKernel : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFELiftingFunctionKernel(const InputParameters & parameters);

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
#endif // DWARFELEPHANTFELIFTINGFUNCTIONKERNEL_H
