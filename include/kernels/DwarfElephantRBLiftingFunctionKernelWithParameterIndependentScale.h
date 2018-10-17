/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBLIFTINGFUNCTIONKERNELWITHPARAMETERINDEPENDENTSCALE_H
#define DWARFELEPHANTRBLIFTINGFUNCTIONKERNELWITHPARAMETERINDEPENDENTSCALE_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantRBKernel.h"
#include "Function.h"

///-------------------------------------------------------------------------
class DwarfElephantRBLiftingFunctionKernelWithParameterIndependentScale;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBLiftingFunctionKernelWithParameterIndependentScale>();

///This Kernel is implements the concept of the lifting function to avoid problems caused by non-homogenous DirichletBC. Note that the scale is fixed for all parameters.
class DwarfElephantRBLiftingFunctionKernelWithParameterIndependentScale : public DwarfElephantRBKernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBLiftingFunctionKernelWithParameterIndependentScale(const InputParameters & parameters);

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
#endif // DWARFELEPHANTRBLIFTINGFUNCTIONKERNELWITHPARAMETERINDEPENDENTSCALE_H
