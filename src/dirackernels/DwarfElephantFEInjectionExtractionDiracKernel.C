#include "DwarfElephantFEInjectionExtractionDiracKernel.h"

 registerMooseObject("DwarfElephantApp", DwarfElephantFEInjectionExtractionDiracKernel);

template <>
InputParameters
validParams<DwarfElephantFEInjectionExtractionDiracKernel>()
{
  InputParameters params = validParams<DiracKernel>();
  params.addRequiredParam<Real>("value", "The value of the point source");
  params.addRequiredParam<std::vector<Real>>("point", "The x,y,z coordinates of the point");
  params.addParam<MooseEnum>(
      "source_type", DwarfElephantFEInjectionExtractionDiracKernel::Type() = "injection", "The source type.");
  params.addParam<Real>("start_time", 0.0, "The starting time.");
  params.addParam<Real>("end_time", 0.0, "The ending time.");
  params.addParam<Real>("fluid_density", 0.0, "The density of the fluid.");
  params.declareControllable("value");
  return params;
}

DwarfElephantFEInjectionExtractionDiracKernel::DwarfElephantFEInjectionExtractionDiracKernel(const InputParameters & parameters)
  : DiracKernel(parameters),
    _value(getParam<Real>("value")),
    _point_param(getParam<std::vector<Real>>("point")),
    _source_type(getParam<MooseEnum>("source_type")),
    _start_time(getParam<Real>("start_time")),
    _end_time(getParam<Real>("end_time")),
    _fluid_density(getParam<Real>("fluid_density"))
{
  _p(0) = _point_param[0];

  if (_point_param.size() > 1)
  {
    _p(1) = _point_param[1];

    if (_point_param.size() > 2)
    {
      _p(2) = _point_param[2];
    }
  }
}

void
DwarfElephantFEInjectionExtractionDiracKernel::addPoints()
{
  addPoint(_p);
}

MooseEnum
DwarfElephantFEInjectionExtractionDiracKernel::Type()
{
  return MooseEnum("injection=1 extraction=2");
}

Real
DwarfElephantFEInjectionExtractionDiracKernel::computeQpResidual()
{
  Real pre_factor = 1.0 / _fluid_density;
  if (_source_type == 1)
    pre_factor *= -1;

  if (_t < _start_time || _t - _dt >= _end_time)
    pre_factor = 0.0;
  else if (_t - _dt < _start_time)
  {
    if (_t <= _end_time)
      pre_factor *= (_t - _start_time) / _dt;
    else
      pre_factor *= (_end_time - _start_time) / _dt;
  }
  else if (_t > _end_time)
    pre_factor *= (_end_time - (_t - _dt)) / _dt;

  return pre_factor  * _value * _test[_i][_qp];
}
