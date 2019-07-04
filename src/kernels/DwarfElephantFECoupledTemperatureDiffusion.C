/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#include "DwarfElephantFECoupledTemperatureDiffusion.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFECoupledTemperatureDiffusion);

template<>
InputParameters validParams<DwarfElephantFECoupledTemperatureDiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addParam<Real>("diffusivity", 1.0, "Diffusivity Coefficient");
  params.addParam<FunctionName>("source", 0.0, "Source Term");
  return params;
}

DwarfElephantFECoupledTemperatureDiffusion::DwarfElephantFECoupledTemperatureDiffusion(const InputParameters & parameters) :
    Diffusion(parameters),
    _diffusivity(getParam<Real>("diffusivity"))
{}

Real
DwarfElephantFECoupledTemperatureDiffusion::computeQpResidual()
{
  return _diffusivity*Diffusion::computeQpResidual();
}

Real
DwarfElephantFECoupledTemperatureDiffusion::computeQpJacobian()
{
  return _diffusivity*Diffusion::computeQpJacobian();
}
