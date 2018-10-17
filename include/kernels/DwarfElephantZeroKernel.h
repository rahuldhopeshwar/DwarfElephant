
///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTZEROKERNEL_H
#define DWARFELEPHANTZEROKERNEL_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "Kernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantZeroKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantZeroKernel>();

class DwarfElephantZeroKernel : public Kernel
{

//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantZeroKernel(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTZEROKERNEL_H
