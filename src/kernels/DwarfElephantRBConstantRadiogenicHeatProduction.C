/**
 * This Kernel is implements a constant radiogenic heat production using the
 * Reduced Basis solution.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBConstantRadiogenicHeatProduction.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBConstantRadiogenicHeatProduction);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBConstantRadiogenicHeatProduction>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  params.addClassDescription("The class implements a constant radiogenic heat"
                              "production.");
  params.addRequiredParam<Real>("radiogenic_heat_production", "Defines the value"
                                "of the radiogenic heat production");
  params.addParam<Real>("norm_value", 1.0, "Defines the normalization value.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBConstantRadiogenicHeatProduction::DwarfElephantRBConstantRadiogenicHeatProduction(const InputParameters & parameters) :
  DwarfElephantRBKernel(parameters),
  _radiogenic_heat_production(getParam<Real>("radiogenic_heat_production")),
  _norm_value(getParam<Real>("norm_value"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantRBConstantRadiogenicHeatProduction::computeQpResidual()
{
  return -(_radiogenic_heat_production/_norm_value) * _test[_i][_qp];
}

Real
DwarfElephantRBConstantRadiogenicHeatProduction::computeQpJacobian()
{
  return 0;
}
