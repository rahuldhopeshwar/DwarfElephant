
///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBONEKERNEL_H
#define DWARFELEPHANTRBONEKERNEL_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBKernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBOneKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBOneKernel>();

///-------------------------------------------------------------------------
class DwarfElephantRBOneKernel : public DwarfElephantRBKernel
{

//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBOneKernel(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBONEKERNEL_H
