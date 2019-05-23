/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBLIFTINGFUNCTIONTIMEKERNELWITHPARAMETERINDEPENDENTSCALE_H
#define DWARFELEPHANTRBLIFTINGFUNCTIONTIMEKERNELWITHPARAMETERINDEPENDENTSCALE_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantRBKernel.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBLiftingFunctionTimeKernelWithParameterIndependentScale;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBLiftingFunctionTimeKernelWithParameterIndependentScale>();

///This Kernel is implements the concept of the lifting function to avoid problems caused by non-homogenous DirichletBC.
class DwarfElephantRBLiftingFunctionTimeKernelWithParameterIndependentScale : public DwarfElephantRBKernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBLiftingFunctionTimeKernelWithParameterIndependentScale(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  Function * _lifting_function;
  Real _scale;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBLIFTINGFUNCTIONTIMEKERNELWITHPARAMETERINDEPENDENTSCALE_H
