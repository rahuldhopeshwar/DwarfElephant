/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
//#include "DwarfElephantFEElectricalConduction.h"
#include "DwarfElephantFEElectricalConduction.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFEElectricalConduction);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEElectricalConduction>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("The class implements an electrical conduction \
                              problem.");
  params.addRequiredParam<Real>("resistivity", "The electrical resistivity of the layer.");
  params.addParam<Real>("norm_value", 1.0, "Defines the normalization value.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEElectricalConduction::DwarfElephantFEElectricalConduction(const InputParameters & parameters) :
  Diffusion(parameters),
  _resistivity(getParam<Real>("resistivity")),
  _norm_value(getParam<Real>("norm_value"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEElectricalConduction::computeQpResidual()
{
  return (_norm_value/_resistivity) * Diffusion::computeQpResidual();
}

Real
DwarfElephantFEElectricalConduction::computeQpJacobian()
{
  return (_norm_value/_resistivity) * Diffusion::computeQpJacobian();
}
