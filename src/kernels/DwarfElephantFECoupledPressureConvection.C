/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#include "DwarfElephantFECoupledPressureConvection.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFECoupledPressureConvection);

template<>
InputParameters validParams<DwarfElephantFECoupledPressureConvection>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("pressure", "Pressure is required for DwarfElephantFECoupledPressureConvection.");
  params.addParam<unsigned>("evaluation_component", "0: Rayleigh number influence is evaulated only in x-direction,"
                                                    "1: Rayleigh number influence is evaulated only in y-direction,"
                                                    "2: Rayleigh number influence is evaulated only in z-direction");
  return params;
}

DwarfElephantFECoupledPressureConvection::DwarfElephantFECoupledPressureConvection(const InputParameters & parameters) :
    Kernel(parameters),
    _grad_p(coupledGradient("pressure")),
    _p(coupledValue("pressure")),
    _second_temp(coupledSecond("pressure")),
    _second_u(second()),
    _second_test(secondTest()),
    _second_phi(secondPhi()),
    _grad_p_var_num(coupled("pressure")),
    _Ra(getMaterialProperty<Real>("rayleigh_material")),
    _component(getParam<unsigned>("evaluation_component"))
{}

Real DwarfElephantFECoupledPressureConvection::computeQpResidual()
{
  return -_test[_i][_qp]*_Ra[_qp]*(_grad_p[_qp]*_grad_u[_qp])
          + _test[_i][_qp]*_Ra[_qp]*_Ra[_qp]*_u[_qp]*_grad_u[_qp](_component);
}

Real DwarfElephantFECoupledPressureConvection::computeQpJacobian()
{
  return -_test[_i][_qp]*_Ra[_qp]*(_grad_p[_qp]*_grad_phi[_j][_qp])
         + _test[_i][_qp]*_Ra[_qp]*_Ra[_qp]*_phi[_j][_qp]*_grad_u[_qp](_component)
         + _test[_i][_qp]*_Ra[_qp]*_Ra[_qp]*_u[_qp]*_grad_phi[_j][_qp](_component);
}

Real DwarfElephantFECoupledPressureConvection::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _grad_p_var_num)
    return -_test[_i][_qp]*_Ra[_qp]*(_grad_phi[_j][_qp]*_grad_u[_qp]);
  else
    return 0;
}
