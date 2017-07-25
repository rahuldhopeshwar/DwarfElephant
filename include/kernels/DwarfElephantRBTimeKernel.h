#ifndef DWARFELEPHANTRBTIMEKERNEL_H
#define DWARFELEPHANTRBTIMEKERNEL_H

#include "DwarfElephantRBKernel.h"

// Forward Declaration
class DwarfElephantRBTimeKernel;

template <>
InputParameters validParams<TimeKernel>();

/**
 * All time kernels should inherit from this class
 *
 */
class DwarfElephantRBTimeKernel : public DwarfElephantRBKernel
{
public:
  DwarfElephantRBTimeKernel(const InputParameters & parameters);

  virtual void computeResidual() override;
};

#endif // DWARFELEPHANTRBTIMEKERNEL_H
