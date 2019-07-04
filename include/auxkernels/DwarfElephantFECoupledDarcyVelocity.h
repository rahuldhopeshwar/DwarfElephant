/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#ifndef DWARFELEPHANTFECOUPLEDDARCYVELOCITY_H
#define DWARFELEPHANTFECOUPLEDDARCYVELOCITY_H

#include "AuxKernel.h"


//Forward Declarations
class DwarfElephantFECoupledDarcyVelocity;

template<>
InputParameters validParams<DwarfElephantFECoupledDarcyVelocity>();

/**
 * Coupled auxiliary value
 */
class DwarfElephantFECoupledDarcyVelocity : public AuxKernel
{
public:
  DwarfElephantFECoupledDarcyVelocity(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const VariableGradient & _grad_pressure;
  const VariableValue & _temp;
  unsigned _component;
  const MaterialProperty<Real> & _Ra;
  Real _temp_ref;
};

#endif //DWARFELEPHANTFECOUPLEDDARCYVELOCITY_H
