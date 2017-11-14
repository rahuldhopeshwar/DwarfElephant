#ifndef DWARFELEPHANTRBINTEGRATEDBC_H
#define DWARFELEPHANTRBINTEGRATEDBC_H

// MOOSE includes
#include "IntegratedBC.h"
#include "BlockRestrictable.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"
#include "MooseVariableScalar.h"

// Forward declarations
//class BlockRestrictable;

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

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOutput();
  virtual void computeJacobianBlock(unsigned int jvar) override;
  void computeJacobianBlockScalar(unsigned int jvar);
  virtual void computeNonlocalJacobian() override {}
  virtual void computeNonlocalOffDiagJacobian(unsigned int /* jvar */) override {}
  virtual void initialSetup() override;

protected:

  /*Methods*/
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /*Attributes*/
  bool _use_displaced;
  bool _matrix_seperation_according_to_subdomains;
  bool _compute_output;

  std::string _simulation_type;

  unsigned int _ID_first_block;
  unsigned int _ID_Aq;
  unsigned int _ID_Mq;
  unsigned int _ID_Fq;
  unsigned int _ID_Oq;

  Real _max_x;
  Real _min_x;
  Real _max_y;
  Real _min_y;
  Real _max_z;
  Real _min_z;
  Real _output_volume;

  DenseVector<Number> _local_out;
  EquationSystems & _es;
};

#endif /* DWARFELEPHANTRBINTEGRATEDBC_H */
