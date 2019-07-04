/* This class was taken from the MOOSE Application beagle written by Powei Huang.
   We transferred it to this package to ensure that all classes are running with
   the same MOOSE version. */

#ifndef DWARFELEPHANTFECOUPLEDPRESSURECONVECTION_H
#define DWARFELEPHANTFECOUPLEDPRESSURECONVECTION_H

#include "Kernel.h"

class DwarfElephantFECoupledPressureConvection;

template<>
InputParameters validParams<DwarfElephantFECoupledPressureConvection>();

class DwarfElephantFECoupledPressureConvection : public Kernel
{
public:

  DwarfElephantFECoupledPressureConvection(const InputParameters & parameters);
  virtual ~DwarfElephantFECoupledPressureConvection() {}

protected:

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  const VariableGradient & _grad_p;
  const VariableValue & _p;
  const VariableSecond & _second_temp;
  const VariableSecond & _second_u;
  const VariableTestSecond & _second_test;
  const VariablePhiSecond & _second_phi;
  unsigned _grad_p_var_num;
  const MaterialProperty<Real> & _Ra;
  const unsigned _component;
};

#endif //DwarfElephantFECoupledPressureConvection_H
