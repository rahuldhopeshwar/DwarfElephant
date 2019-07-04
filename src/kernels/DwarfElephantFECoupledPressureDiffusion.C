/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#include "DwarfElephantFECoupledPressureDiffusion.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFECoupledPressureDiffusion);

template<>
InputParameters validParams<DwarfElephantFECoupledPressureDiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addCoupledVar("temperature","temperature is required for DwarfElephantFECoupledPressureDiffusion.");
  params.addParam<unsigned>("evaluation_component", "0: Kernel is evaulated only in x-direction,"
                                                    "1: Kernel is evaulated only in y-direction,"
                                                    "2: Kernel is evaulated only in z-direction");
  return params;
}

DwarfElephantFECoupledPressureDiffusion::DwarfElephantFECoupledPressureDiffusion(const InputParameters & parameters) :
    Diffusion(parameters),
    _temp(coupledValue("temperature")),
    _temp_var_num(coupled("temperature")),
    _grad_temp(coupledGradient("temperature")),
    _Ra(getMaterialProperty<Real>("rayleigh_material")),
    _component(getParam<unsigned>("evaluation_component"))
{}

Real
DwarfElephantFECoupledPressureDiffusion::computeQpResidual()
{
  return Diffusion::computeQpResidual() - _Ra[_qp]*_grad_test[_i][_qp](_component)*_temp[_qp];
}

Real
DwarfElephantFECoupledPressureDiffusion::computeQpJacobian()
{
  return Diffusion::computeQpJacobian();
}


Real
DwarfElephantFECoupledPressureDiffusion::computeQpOffDiagJacobian(unsigned jvar)
{
  if(jvar == _temp_var_num)
    return  -_Ra[_qp]*_grad_test[_i][_qp](_component)*_phi[_j][_qp];
  else
    return 0;
}
