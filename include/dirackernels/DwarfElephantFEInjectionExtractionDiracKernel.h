#ifndef DWARFELEPHANTFEINJECTIONEXTRACTIONDIRACKERNEL_H
#define DWARFELEPHANTFEINJECTIONEXTRACTIONDIRACKERNEL_H

// Moose Includes
#include "DiracKernel.h"

// Forward Declarations
class DwarfElephantFEInjectionExtractionDiracKernel;

template <>
InputParameters validParams<DwarfElephantFEInjectionExtractionDiracKernel>();

///Analoge to the class ConstantPointSource. It enables the use of RB for this DiracKernel. Note, that the value of the point source is fixed for all parameters.
class DwarfElephantFEInjectionExtractionDiracKernel : public DiracKernel
{
public:
  DwarfElephantFEInjectionExtractionDiracKernel(const InputParameters & parameters);

  virtual void addPoints() override;
  static MooseEnum Type();

protected:
  virtual Real computeQpResidual() override;

  const Real & _value;
  std::vector<Real> _point_param;
  Point _p;
  MooseEnum _source_type;
  Real _start_time;
  Real _end_time;
  Real _fluid_density;
};

#endif // DWARFELEPHANTFEINJECTIONEXTRACTIONDIRACKERNEL_H
