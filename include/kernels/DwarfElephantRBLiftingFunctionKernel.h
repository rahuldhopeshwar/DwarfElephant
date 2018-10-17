/**
 * This Kernel is implements the concept of the lifting function to avoid
 * problems caused by non-homogenous DirichletBC.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBLIFTINGFUNCTIONKERNEL_H
#define DWARFELEPHANTRBLIFTINGFUNCTIONKERNEL_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantRBKernel.h"
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBLiftingFunctionKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBLiftingFunctionKernel>();

///This Kernel is implements the concept of the lifting function to avoid problems caused by non-homogenous DirichletBC.
class DwarfElephantRBLiftingFunctionKernel : public DwarfElephantRBKernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBLiftingFunctionKernel(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  Function * _lifting_function;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBLIFTINGFUNCTIONKERNEL_H
