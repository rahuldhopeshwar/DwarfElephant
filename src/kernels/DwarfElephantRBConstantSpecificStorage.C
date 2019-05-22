/**
 * This Kernel is implements a constant radiogenic heat production using the
 * Reduced Basis solution.
 */

 //---------------------------------INCLUDES-------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBConstantSpecificStorage.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBConstantSpecificStorage);

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBConstantSpecificStorage>()
{
  InputParameters params = validParams<DwarfElephantRBTimeDerivative>();
  params.addClassDescription("The class implements a constant radiogenic heat"
                              "production.");
  params.addRequiredParam<Real>("specific_storage", "Defines the value"
                                "of the specific storage");
  params.addParam<Real>("norm_value", 1.0, "Defines the normalization value.");
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBConstantSpecificStorage::DwarfElephantRBConstantSpecificStorage(const InputParameters & parameters) :
  DwarfElephantRBTimeDerivative(parameters),
  _specific_storage(getParam<Real>("specific_storage")),
  _norm_value(getParam<Real>("norm_value"))
{
}

//----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
Real
DwarfElephantRBConstantSpecificStorage::computeQpResidual()
{
  return 0;
}

Real
DwarfElephantRBConstantSpecificStorage::computeQpJacobian()
{
  return (_specific_storage/_norm_value) * _phi[_j][_qp] * _test[_i][_qp];
}
