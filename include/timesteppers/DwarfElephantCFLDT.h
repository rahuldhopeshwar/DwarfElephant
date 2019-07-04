/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#ifndef DWARFELEPHANTCFLDT_H
#define DWARFELEPHANTCFLDT_H

#include "TimeStepper.h"
#include "PostprocessorInterface.h"

class DwarfElephantCFLDT;

template <>
InputParameters validParams<DwarfElephantCFLDT>();

class DwarfElephantCFLDT : public TimeStepper, public PostprocessorInterface
{
public:
  DwarfElephantCFLDT(const InputParameters & parameters);

protected:
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;

  const PostprocessorValue & _pps_value;
  const Real & _max_Ra;

  bool _has_initial_dt;
  Real _initial_dt;

  /// Multiplier applied to the postprocessor value
  const Real & _cfl_num;
  const Real & _activate_time;

  /// Factor added to the postprocessor value
  const Real & _factor;
};

#endif /* DwarfElephantCFLDT_H */
