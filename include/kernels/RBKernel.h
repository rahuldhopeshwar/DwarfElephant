/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBKernel inherits from the Kernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * whole stiffness matrix is needed and not only the diagonal entries.
 */

///-------------------------------------------------------------------------
#ifndef RBKERNEL_H
#define RBKERNEL_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Kernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class RBKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<RBKernel>();

///-------------------------------------------------------------------------
class RBKernel : public Kernel
{

//----------------------------------PUBLIC----------------------------------
public:
  RBKernel(const InputParameters & parameters);

 /* Methods */
  virtual void computeResidual() override;
  virtual void computeJacobian() override;

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

//  std::ofstream _file_output;
//  _file_output.open("test_file.txt");

  friend class RBOutput;
  };

///-------------------------------------------------------------------------
#endif //RBKERNEL_H
