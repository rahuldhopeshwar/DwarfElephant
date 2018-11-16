/**
 * This Kernel is implements a constant radiogenic heat production using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEConstantRadiogenicHeatProduction.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFEConstantRadiogenicHeatProduction);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEConstantRadiogenicHeatProduction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The class implements a constant radiogenic heat"
                              "production.");
  params.addRequiredParam<Real>("radiogenic_heat_production", "Defines the value"
                                "of the radiogenic heat production");
  params.addParam<Real>("norm_value", 1.0, "Defines the normalization value.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEConstantRadiogenicHeatProduction::DwarfElephantFEConstantRadiogenicHeatProduction(const InputParameters & parameters) :
  Kernel(parameters),
  _radiogenic_heat_production(getParam<Real>("radiogenic_heat_production")),
  _norm_value(getParam<Real>("norm_value"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEConstantRadiogenicHeatProduction::computeQpResidual()
{
  return -(_radiogenic_heat_production/_norm_value) * _test[_i][_qp];
}

Real
DwarfElephantFEConstantRadiogenicHeatProduction::computeQpJacobian()
{
  return 0;
}
