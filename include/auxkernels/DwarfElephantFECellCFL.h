/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#ifndef DWARFELEPHANTFECELLCFL_H
#define DWARFELEPHANTFECELLCFL_H

// MOOSE includes
#include "AuxKernel.h"

// Forward declarations
class DwarfElephantFECellCFL;

template <>
InputParameters validParams<DwarfElephantFECellCFL>();

class DwarfElephantFECellCFL : public AuxKernel
{
public:
  DwarfElephantFECellCFL(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  FEProblemBase & _feproblem;

private:
  const VariableValue & _velocity_x;
  const VariableValue & _velocity_y;
  const VariableValue & _velocity_z;

  const MaterialProperty<Real> & _scale;
};

#endif // DwarfElephantFECellCFL_H
