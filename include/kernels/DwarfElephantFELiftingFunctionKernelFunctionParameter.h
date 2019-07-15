/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFELIFTINGFUNCTIONKERNELFUNCTIONPARAMETER_H
#define DWARFELEPHANTFELIFTINGFUNCTIONKERNELFUNCTIONPARAMETER_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFELiftingFunctionKernelFunctionParameter;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFELiftingFunctionKernelFunctionParameter>();

///This Kernel is implements the concept of the lifting function to avoid problems caused by non-homogenous DirichletBC.
class DwarfElephantFELiftingFunctionKernelFunctionParameter : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFELiftingFunctionKernelFunctionParameter(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  const Function * _lifting_function;
  Real _scale;
  Real _norm_value;

  const Function & _func;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFELIFTINGFUNCTIONKERNEL_H
