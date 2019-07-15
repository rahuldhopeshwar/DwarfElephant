/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBLIFTINGFUNCTIONKERNELFUNCTIONPARAMETER_H
#define DWARFELEPHANTRBLIFTINGFUNCTIONKERNELFUNCTIONPARAMETER_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantRBKernel.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBLiftingFunctionKernelFunctionParameter;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBLiftingFunctionKernelFunctionParameter>();

///This Kernel is implements the concept of the lifting function to avoid problems caused by non-homogenous DirichletBC.
class DwarfElephantRBLiftingFunctionKernelFunctionParameter : public DwarfElephantRBKernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBLiftingFunctionKernelFunctionParameter(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  const Function * _lifting_function;
  const Function & _func;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBLIFTINGFUNCTIONKERNELFUNCTIONPARAMETER_H
