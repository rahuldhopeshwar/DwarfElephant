/**
 * This Kernel is implements a constant radiogenic heat production using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantFEConstantSpecificStorage.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFEConstantSpecificStorage);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEConstantSpecificStorage>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addClassDescription("The class implements a constant radiogenic heat"
                              "production.");
  params.addRequiredParam<Real>("specific_storage", "Defines the value"
                                "of the specific storage");
  params.addParam<Real>("norm_value", 1.0, "Defines the normalization value.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantFEConstantSpecificStorage::DwarfElephantFEConstantSpecificStorage(const InputParameters & parameters) :
  TimeDerivative(parameters),
  _specific_storage(getParam<Real>("specific_storage")),
  _norm_value(getParam<Real>("norm_value"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantFEConstantSpecificStorage::computeQpResidual()
{
  return (_specific_storage/_norm_value) * _u_dot[_qp] * _test[_i][_qp];
}

Real
DwarfElephantFEConstantSpecificStorage::computeQpJacobian()
{
  return (_specific_storage/_norm_value) * _du_dot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}
