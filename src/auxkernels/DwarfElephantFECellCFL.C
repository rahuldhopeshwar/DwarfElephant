/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#include "DwarfElephantFECellCFL.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFECellCFL);

template <>
InputParameters
validParams<DwarfElephantFECellCFL>()
{
  MooseEnum component("x=0 y=1 z=2");
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("velocity_x", "The velocity in x component");
  params.addRequiredCoupledVar("velocity_y", "The velocity in y component");
  params.addRequiredCoupledVar("velocity_z", "The velocity in z component");
  return params;
}

DwarfElephantFECellCFL::DwarfElephantFECellCFL(const InputParameters & parameters)
  : AuxKernel(parameters),
    _feproblem(dynamic_cast<FEProblemBase &>(_subproblem)),
    _velocity_x(coupledValue("velocity_x")),
    _velocity_y(coupledValue("velocity_y")),
    _velocity_z(coupledValue("velocity_z")),
    _scale(getMaterialProperty<Real>("rayleigh_material"))
{
}

Real
DwarfElephantFECellCFL::computeValue()
{
  return _scale[_qp]*(std::abs(_velocity_x[_qp]) + std::abs(_velocity_y[_qp]) + std::abs(_velocity_z[_qp]))*_feproblem.dt()/_current_elem->hmin();
}
