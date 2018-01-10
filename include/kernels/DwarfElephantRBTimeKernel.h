/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBTimeKernel inherits from the RBKernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * mass matrix is needed separated in its subdomain contributions. In
 * addition it overwrites the function computeResidual.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBTIMEKERNEL_H
#define DWARFELEPHANTRBTIMEKERNEL_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBKernel.h"

///-------------------------------------------------------------------------
// Forward Declaration
class DwarfElephantRBTimeKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template <>
InputParameters validParams<TimeKernel>();

///-------------------------------------------------------------------------
/**
 * All time kernels should inherit from this class
 *
 */
class DwarfElephantRBTimeKernel : public DwarfElephantRBKernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBTimeKernel(const InputParameters & parameters);

  /*Methods*/
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBTIMEKERNEL_H
