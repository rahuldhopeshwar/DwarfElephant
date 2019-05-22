#include "DwarfElephantRBPointSource.h"

 registerMooseObject("DwarfElephantApp", DwarfElephantRBPointSource);

template <>
InputParameters
validParams<DwarfElephantRBPointSource>()
{
  InputParameters params = validParams<DwarfElephantRBDiracKernel>();
  params.addRequiredParam<std::vector<Real>>("point", "The x,y,z coordinates of the point");
  params.addParam<bool>("sink", false, "Defines that the source term is treated as a sink term.");
  return params;
}

DwarfElephantRBPointSource::DwarfElephantRBPointSource(const InputParameters & parameters)
  : DwarfElephantRBDiracKernel(parameters),
    _point_param(getParam<std::vector<Real>>("point")),
    _sink(getParam<bool>("sink"))
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
DwarfElephantRBPointSource::addPoints()
{
  addPoint(_p);
}

Real
DwarfElephantRBPointSource::computeQpResidual()
{
  Real scale = 1.0;
  if(_sink)
    scale *=-1.0;
  //  This is negative because it's a forcing function that has been brought over to the left side
  return scale * _test[_i][_qp];
}
