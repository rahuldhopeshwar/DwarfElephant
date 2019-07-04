/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#include "DwarfElephantCFLDT.h"
#include "FEProblem.h"

registerMooseObject("DwarfElephantApp", DwarfElephantCFLDT);

template <>
InputParameters
validParams<DwarfElephantCFLDT>()
{
  InputParameters params = validParams<TimeStepper>();
  params.addRequiredParam<PostprocessorName>("postprocessor",
                                             "The name of the postprocessor that computes the dt");
  params.addRequiredParam<Real>("max_Ra",
                                "Maximum sqrt(Ra)");
  params.addParam<Real>("dt", "Initial value of dt");
  params.addParam<Real>("cfl", 1, "Desired CFL value.");
  params.addParam<Real>("activate_time", "Desired time to start DwarfElephantCFLDT");
  params.addParam<Real>("factor", 0, "Add a factor to the supplied postprocessor value.");
  return params;
}

DwarfElephantCFLDT::DwarfElephantCFLDT(const InputParameters & parameters)
  : TimeStepper(parameters),
    PostprocessorInterface(this),
    _pps_value(getPostprocessorValue("postprocessor")),
    _max_Ra(getParam<Real>("max_Ra")),
    _has_initial_dt(isParamValid("dt")),
    _initial_dt(_has_initial_dt ? getParam<Real>("dt") : 0.),
    _cfl_num(getParam<Real>("cfl")),
    _activate_time(getParam<Real>("activate_time")),
    _factor(getParam<Real>("factor"))
{
}

Real
DwarfElephantCFLDT::computeInitialDT()
{
  if (_has_initial_dt)
    return _initial_dt;
  else
    return computeDT();
}

Real
DwarfElephantCFLDT::computeDT()
{
  if (_time < _activate_time)
    return _initial_dt;
  else
    return _cfl_num / _max_Ra * _pps_value + _factor;
}
