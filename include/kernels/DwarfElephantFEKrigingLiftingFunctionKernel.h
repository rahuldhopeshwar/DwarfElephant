/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFEKRIGINGLIFTINGFUNCTIONKERNEL_H
#define DWARFELEPHANTFEKRIGINGLIFTINGFUNCTIONKERNEL_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEKrigingLiftingFunctionKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEKrigingLiftingFunctionKernel>();

///This Kernel is implements the concept of the lifting function to avoid problems caused by non-homogenous DirichletBC.
class DwarfElephantFEKrigingLiftingFunctionKernel : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEKrigingLiftingFunctionKernel(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  Real _range;
  const Function & _lifting_function_1;
  const Function & _lifting_function_2;
  Real _scale;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFEKRIGINGLIFTINGFUNCTIONKERNEL_H
