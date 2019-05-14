#ifndef DWARFELEPHANTRBITERATIONADAPTIVEDT_H
#define DWARFELEPHANTRBITERATIONADAPTIVEDT_H

#include "DwarfElephantRBExecutioner.h"
#include "DwarfElephantInitializeRBSystemTransient.h"

class DwarfElephantRBIterationAdaptiveDT;
class DwarfElephantInitializeRBSystemTransient;

template <>
InputParameters validParams<DwarfElephantRBIterationAdaptiveDT>();

/**
 * All Distributions should inherit from this class
 */
class DwarfElephantRBIterationAdaptiveDT : public DwarfElephantRBExecutioner
{
public:
  DwarfElephantRBIterationAdaptiveDT(const InputParameters & parameters);
  /**
   * Compute the timestep
   */
  void execute() override;
  virtual Real computeDT() override;

  unsigned int _n_l_iter;

protected:
  /// Adapt the timestep to maintain this non-linear iteration count...
  int _optimal_iterations;
  /// ...plus/minus this value.
  int _iteration_window;
  /// use _optimal_iterations and _iteration_window multiplied with this factor for linear iterations
  const int _linear_iteration_ratio;
  /// grow the timestep by this factor
  const Real & _growth_factor;
  /// cut the timestep by by this factor
  const Real & _cutback_factor;
  bool & _cutback_occurred;
  Real dt;
};

#endif /* DWARFELEPHANTRBITERATIONADAPTIVEDT_H */
