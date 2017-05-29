#ifndef RBINTEGRATEDBC_H
#define RBINTEGRATEDBC_H

// MOOSE includes
#include "IntegratedBC.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystem.h"

// Forward declarations
class DwarfElephantInitializeRBSystem;
class RBIntegratedBC;

template<>
InputParameters validParams<RBIntegratedBC>();

class RBIntegratedBC :
  public IntegratedBC
{
public:

  /*Methods*/
  RBIntegratedBC(const InputParameters & parameters);

  virtual ~RBIntegratedBC();

  virtual void computeResidual();
  virtual void computeJacobian();
  virtual void computeJacobianBlock(unsigned int jvar);
  void computeJacobianBlockScalar(unsigned int jvar);
  virtual void computeNonlocalJacobian() {}
  virtual void computeNonlocalOffDiagJacobian(unsigned int /* jvar */) {}

protected:

  /*Methods*/
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /*Attributes*/
  bool _use_displaced;

  Real _output_volume;
  DenseVector<Number> _local_out;

  EquationSystems & _es;
  const DwarfElephantInitializeRBSystem & _initialize_rb_system;
};

#endif /* RBINTEGRATEDBC_H */
