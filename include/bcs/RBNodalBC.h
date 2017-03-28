#ifndef RBNODALBC_H
#define RBNODALBC_H

#include "NodalBC.h"
#include "MooseMesh.h"
#include "NonlinearSystemBase.h"

#include "DwarfElephantInitializeRBSystem.h"
#include "CacheStiffnessMatrix.h"

// Forward declarations
class MooseMesh;
class NonlinearSystemBase;
class DwarfElephantInitializeRBSystem;
class RBNodalBC;

template<>
InputParameters validParams<RBNodalBC>();

class RBNodalBC :
  public NodalBC
{
public:
  RBNodalBC(const InputParameters & parameters);

  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const DwarfElephantInitializeRBSystem & _initialize_rb_system;
};

#endif /* RBNODALBC_H */
