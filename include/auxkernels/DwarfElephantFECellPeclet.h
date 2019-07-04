/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#ifndef DWARFELEPHANTFECELLPECLET_H
#define DWARFELEPHANTFECELLPECLET_H

// MOOSE includes
#include "AuxKernel.h"

// Forward declarations
class DwarfElephantFECellPeclet;

template <>
InputParameters validParams<DwarfElephantFECellPeclet>();

/**
 * Extract a component from the gradient of a variable
 */
class DwarfElephantFECellPeclet : public AuxKernel
{
public:
  DwarfElephantFECellPeclet(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  const VariableValue & _velocity_x;
  const VariableValue & _velocity_y;
  const VariableValue & _velocity_z;

  const MaterialProperty<Real> & _scale;
};

#endif // DwarfElephantFECellPeclet_H
