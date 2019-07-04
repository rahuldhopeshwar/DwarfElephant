/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#include "DwarfElephantFECoupledDarcyVelocity.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFECoupledDarcyVelocity);

template<>
InputParameters validParams<DwarfElephantFECoupledDarcyVelocity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("pressure", "pressure");
  params.addCoupledVar("temperature", "temperature");
  params.addRequiredParam<unsigned>("evaluation_component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");

  return params;
}

DwarfElephantFECoupledDarcyVelocity::DwarfElephantFECoupledDarcyVelocity(const InputParameters & parameters) :
    AuxKernel(parameters),
    _grad_pressure(coupledGradient("pressure")),
    _temp(coupledValue("temperature")),
    _component(getParam<unsigned>("evaluation_component")),
    _Ra(getMaterialProperty<Real>("rayleigh_material"))
{}

Real
DwarfElephantFECoupledDarcyVelocity::computeValue()
{
  return -_grad_pressure[_qp](_component) + _temp[_qp];
}
