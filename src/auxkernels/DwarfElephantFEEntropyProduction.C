/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   It has been modified to be applicable for general geometries.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#include "DwarfElephantFEEntropyProduction.h"

registerMooseObject("DwarfElephantApp", DwarfElephantFEEntropyProduction);

template<>
InputParameters validParams<DwarfElephantFEEntropyProduction>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredParam<FunctionName>("T_top", "A function that describes the temperature at the top of the model.");
  params.addRequiredParam<FunctionName>("T_bottom", "A function that describes the temperature at the bottom of the model.");
  params.addRequiredParam<FunctionName>("z_top", "A function that describes the z-coordinates at the top of the model.");
  params.addRequiredParam<FunctionName>("z_bottom", "A function that describes the z-coordinates at the bottom of the model.");
  params.addParam<Real>("gravity_acceleration", 9.81, "The gravity acceleration (default value: 9.81 m/s^2).");
  params.addParam<Real>("alpha", "Thermal expansion coefficient (T^-1)");
  params.addParam<Real>("cf", "Heat capacity of the fluid");
  params.addParam<Real>("thermal_conductivity", 0.0, "Thermal conductivity of the unit.");
  params.addParam<bool>("entropy_generation_number", true, "If true the entropy generation number is calculated else the entropy production.");
  params.addRequiredCoupledVar("temp", "Entropy Production AuxKernel requires temperature");
  params.addRequiredCoupledVar("velocity_x", "Entropy Production AuxKernel requires velocity_x");
  params.addCoupledVar("velocity_y", "Entropy Production AuxKernel: velocity_y");
  params.addCoupledVar("velocity_z", "Entropy Production AuxKernel: velocity_z");
  return params;
}

DwarfElephantFEEntropyProduction::DwarfElephantFEEntropyProduction(const InputParameters & parameters) :
    AuxKernel(parameters),
    _T_top(getFunction("T_top")),
    _T_bottom(getFunction("T_bottom")),
    _z_top(getFunction("z_top")),
    _z_bottom(getFunction("z_bottom")),
    _gravity_acceleration(getParam<Real>("gravity_acceleration")),
    _alpha(getParam<Real>("alpha")),
    _cf(getParam<Real>("cf")),
    _lambda(getParam<Real>("thermal_conductivity")),
    _entropy_generation_number(getParam<bool>("entropy_generation_number")),
    _grad_temp(coupledGradient("temp")),
    _temp(coupledValue("temp")),
    _vel_x(coupledValue("velocity_x")),
    _vel_y(coupledValue("velocity_y")),
    _vel_z(coupledValue("velocity_z"))
{
  if(!_entropy_generation_number && _lambda == 0)
    mooseError("You have to define a thermal conducitivity.");
}

Real
DwarfElephantFEEntropyProduction::computeValue()
{
  Real _T_bar = (_T_top.value(_c_fe_problem.time(),_q_point[_qp])+
                 _T_bottom.value(_c_fe_problem.time(),_q_point[_qp]))/2;
  Real _delta_T = _T_top.value(_c_fe_problem.time(),_q_point[_qp])+
                  _T_bottom.value(_c_fe_problem.time(),_q_point[_qp]);

  Real _distance = abs(_z_top.value(_c_fe_problem.time(),_q_point[_qp])-
                   _z_bottom.value(_c_fe_problem.time(),_q_point[_qp]));

  Real _q_square = _vel_x[_qp]*_vel_x[_qp] + _vel_y[_qp]*_vel_y[_qp] +
                   _vel_z[_qp]*_vel_z[_qp];

  Real _entropy_therm_gen_num = _grad_temp[_qp]*_grad_temp[_qp];
  Real _entropy_visc_gen_num = _alpha*_distance*_T_bar*_gravity_acceleration/
                               (_delta_T*_cf)*(_q_square);

  Real _entropy_factor = (_lambda*_delta_T*_delta_T)/
                         (_distance*_distance*_T_bar*_T_bar);

  if(_entropy_generation_number)
    return _entropy_therm_gen_num + _entropy_visc_gen_num;
  else
    return _entropy_factor*(_entropy_therm_gen_num + _entropy_visc_gen_num);
}
