#include "DwarfElephantRBConstantPointSource.h"

template <>
InputParameters
validParams<DwarfElephantRBConstantPointSource>()
{
  InputParameters params = validParams<DwarfElephantRBDiracKernel>();
  params.addRequiredParam<Real>("value", "The value of the point source");
  params.addRequiredParam<std::vector<Real>>("point", "The x,y,z coordinates of the point");
  params.declareControllable("value");
  return params;
}

DwarfElephantRBConstantPointSource::DwarfElephantRBConstantPointSource(const InputParameters & parameters)
  : DwarfElephantRBDiracKernel(parameters),
    _value(getParam<Real>("value")),
    _point_param(getParam<std::vector<Real>>("point"))
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
DwarfElephantRBConstantPointSource::addPoints()
{
  addPoint(_p);
}

Real
DwarfElephantRBConstantPointSource::computeQpResidual()
{
  //  This is negative because it's a forcing function that has been brought over to the left side
  return _test[_i][_qp] * _value;
}
