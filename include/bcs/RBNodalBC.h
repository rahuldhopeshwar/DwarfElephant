#ifndef RBNODALBC_H
#define RBNODALBC_H

#include "libmesh/equation_systems.h"

#include "NodalBC.h"
#include "MooseMesh.h"
#include "NonlinearSystemBase.h"

#include "DwarfElephantInitializeRBSystem.h"
#include "CacheStiffnessMatrix.h"

// Forward declarations
// libMesh includes
namespace libMesh
{
  class EquationSystems;
}

// MOOSE includes
class MooseMesh;
class NonlinearSystemBase;
class CacheStiffnessMatrix;

// MOOSE includes (DwarfElephant package)
class DwarfElephantInitializeRBSystem;
class RBNodalBC;

template<>
InputParameters validParams<RBNodalBC>();

class RBNodalBC :
  public NodalBC
{
public:
  RBNodalBC(const InputParameters & parameters);

  virtual void computeResidual(NumericVector<Number> & residual) override;
  virtual void computeJacobian() override;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  const DwarfElephantInitializeRBSystem & _initialize_rb_system;

  Function * _function;
  CacheStiffnessMatrix * _cache_stiffness_matrix;
};

#endif /* RBNODALBC_H */
