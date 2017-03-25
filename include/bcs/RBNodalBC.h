#ifndef RBNODALBC_H
#define RBNODALBC_H

#include "NodalBC.h"
#include "MooseMesh.h"
#include "BlockRestrictable.h"

#include "DwarfElephantInitializeRBSystem.h"

// Forward declarations
class MooseMesh;
class DwarfElephantInitializeRBSystem;
class RBNodalBC;

template<>
InputParameters validParams<RBNodalBC>();

class RBNodalBC :
  public NodalBC,
  public BlockRestrictable
{
public:
  RBNodalBC(const InputParameters & parameters);

  virtual void computeJacobian();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
  const DwarfElephantInitializeRBSystem & _initialize_rb_system;
};

#endif /* RBNODALBC_H */
