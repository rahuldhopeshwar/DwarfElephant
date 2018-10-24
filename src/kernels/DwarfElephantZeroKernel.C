 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantZeroKernel.h"

registerMooseObject("DwarfElephantApp", DwarfElephantZeroKernel);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantZeroKernel>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantZeroKernel::DwarfElephantZeroKernel(const InputParameters & parameters) :
  Kernel(parameters)
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantZeroKernel::computeQpResidual()
{
  return 0;
}

Real
DwarfElephantZeroKernel::computeQpJacobian()
{
  return 0;
}
