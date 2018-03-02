#ifndef DWARFELEPHANTRBCONSTANTPOINTSOURCE_H
#define DWARFELEPHANTRBCONSTANTPOINTSOURCE_H

// Moose Includes
#include "DwarfElephantRBDiracKernel.h"

// Forward Declarations
class DwarfElephantRBConstantPointSource;

template <>
InputParameters validParams<DwarfElephantRBConstantPointSource>();

class DwarfElephantRBConstantPointSource : public DwarfElephantRBDiracKernel
{
public:
  DwarfElephantRBConstantPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;

  const Real & _value;
  std::vector<Real> _point_param;
  Point _p;
};

#endif // DWARFELEPHANTRBCONSTANTPOINTSOURCE_H
