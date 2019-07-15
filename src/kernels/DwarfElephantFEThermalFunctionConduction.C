/**
 * This Kernel is implements a thermal conduction problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEThermalFunctionConduction.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFEThermalFunctionConduction);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEThermalFunctionConduction>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription("The class implements a thermal conduction \
                              problem.");
  params.addRequiredParam<Real>("thermal_conductivity", "Defines the value of the thermal conductivity");
  params.addParam<Real>("norm_value", 1.0, "Defines the normalization value.");
  params.addRequiredParam<FunctionName>("function", "Name of the function that describes the parameter dependence.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEThermalFunctionConduction::DwarfElephantFEThermalFunctionConduction(const InputParameters & parameters) :
  Diffusion(parameters),
  _lambda(getParam<Real>("thermal_conductivity")),
  _norm_value(getParam<Real>("norm_value")),
  _func(getFunction("function"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEThermalFunctionConduction::computeQpResidual()
{
  Real dependency = _func.value(_fe_problem.time(),_q_point[_qp]);
  return (_lambda/_norm_value) * dependency * Diffusion::computeQpResidual();
}

Real
DwarfElephantFEThermalFunctionConduction::computeQpJacobian()
{
  Real dependency = _func.value(_fe_problem.time(),_q_point[_qp]);
  return (_lambda/_norm_value) * dependency * Diffusion::computeQpJacobian();
}
