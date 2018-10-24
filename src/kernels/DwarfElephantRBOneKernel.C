 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBOneKernel.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBOneKernel);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBOneKernel>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBOneKernel::DwarfElephantRBOneKernel(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters)
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantRBOneKernel::computeQpResidual()
{
  return -_test[_i][_qp];
}

Real
DwarfElephantRBOneKernel::computeQpJacobian()
{
  return 0;
}
