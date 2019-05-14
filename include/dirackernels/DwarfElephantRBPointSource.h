#ifndef DWARFELEPHANTRBPOINTSOURCE_H
#define DWARFELEPHANTRBPOINTSOURCE_H

// Moose Includes
#include "DwarfElephantRBDiracKernel.h"

// Forward Declarations
class DwarfElephantRBPointSource;

template <>
InputParameters validParams<DwarfElephantRBPointSource>();

///Analoge to the class ConstantPointSource. It enables the use of RB for this DiracKernel. Note, that the value of the point source is fixed for all parameters.
class DwarfElephantRBPointSource : public DwarfElephantRBDiracKernel
{
public:
  DwarfElephantRBPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;

  std::vector<Real> _point_param;
  Point _p;
};

#endif // DWARFELEPHANTRBPOINTSOURCE_H
