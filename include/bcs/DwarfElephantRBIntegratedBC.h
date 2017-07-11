#ifndef DWARFELEPHANTRBINTEGRATEDBC_H
#define DWARFELEPHANTRBINTEGRATEDBC_H

// MOOSE includes
#include "IntegratedBC.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"

// Forward declarations
class DwarfElephantInitializeRBSystemSteadyState;
class DwarfElephantRBIntegratedBC;

template<>
InputParameters validParams<DwarfElephantRBIntegratedBC>();

class DwarfElephantRBIntegratedBC :
  public IntegratedBC
{
public:

  /*Methods*/
  DwarfElephantRBIntegratedBC(const InputParameters & parameters);

  virtual ~DwarfElephantRBIntegratedBC();

  virtual void computeResidual();
  virtual void computeJacobian();
  virtual void computeJacobianBlock(unsigned int jvar);
  void computeJacobianBlockScalar(unsigned int jvar);
  virtual void computeNonlocalJacobian() {}
  virtual void computeNonlocalOffDiagJacobian(unsigned int /* jvar */) {}
  virtual void initialSetup() override;

protected:

  /*Methods*/
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /*Attributes*/
  bool _use_displaced;

  std::string _simulation_type;

  unsigned int _ID_Aq;
  unsigned int _ID_Fq;

  EquationSystems & _es;
};

#endif /* DWARFELEPHANTRBINTEGRATEDBC_H */
