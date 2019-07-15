/* This class was taken from the MOOSE Application beagle written by Powei Huang.
    It has been modified to be applicable for general geometries.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#ifndef DWARFELEPHANTFEENTROPYPRODUCTION_H
#define DWARFELEPHANTFEENTROPYPRODUCTION_H

#include "AuxKernel.h"
#include "Function.h"

// Forward Declarations
class DwarfElephantFEEntropyProduction;

template<>
InputParameters validParams<DwarfElephantFEEntropyProduction>();

class DwarfElephantFEEntropyProduction : public AuxKernel
{
public:
  DwarfElephantFEEntropyProduction(const InputParameters & parameters);

  virtual ~DwarfElephantFEEntropyProduction() {}

protected:
  virtual Real computeValue() override;

  const Function & _T_top;
  const Function & _T_bottom;
  const Function & _z_top;
  const Function & _z_bottom;

  Real _gravity_acceleration;
  Real _alpha;
  Real _cf;
  Real _lambda;

  bool _entropy_generation_number;
  bool _thermal_part;
  bool _visc_part;

  const VariableGradient & _grad_temp;
  const VariableValue & _temp;
  const VariableValue & _vel_x;
  const VariableValue & _vel_y;
  const VariableValue & _vel_z;
};

#endif // DwarfElephantFEEntropyProduction_H
