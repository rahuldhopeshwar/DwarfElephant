#ifndef DWARFELEPHANTCOMPUTEEIMINNERPRODUCTMATRIXSTEADYSTATE_H
#define DWARFELEPHANTCOMPUTEEIMINNERPRODUCTMATRIXSTEADYSTATE_H

#include "ElementUserObject.h"
#include "MooseVariableInterface.h"

#include "libmesh/sparse_matrix.h"

#include "DwarfElephantInitializeRBSystemSteadyState.h"

namespace libMesh
{
  template<typename T> class SparseMatrix;
}

class DwarfElephantComputeEIMInnerProductMatrixSteadyState;

template <>
InputParameters validParams<DwarfElephantComputeEIMInnerProductMatrixSteadyState> ();

class DwarfElephantComputeEIMInnerProductMatrixSteadyState : public ElementUserObject,
                            public MooseVariableInterface<Real>
{
public:
  DwarfElephantComputeEIMInnerProductMatrixSteadyState(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;
  virtual void processRBParameters();

  virtual Real getValue();
  DenseMatrix <Number> _local_ke;

protected:
  virtual Real computeIntegral(unsigned int _i, unsigned int _j);
  
  bool _use_displaced;
  bool _skip_matrix_assembly_in_rb_system;
  bool _skip_vector_assembly_in_rb_system;
  bool _offline_stage;
  bool _compliant;
  bool _deterministic_training_RB;
  bool _quiet_mode_RB;
  bool _normalize_RB_bound_in_greedy;

  unsigned int _n_training_samples_RB;
  unsigned int _training_parameters_random_seed_RB;
  unsigned int _N_max_RB;

  std::string _system_name;
//    std::string _parameters_filename;     //only required if one wants to read the data over the GetPot class from libMesh directly
  std::vector<std::string> _continuous_parameters_RB;
  std::vector<std::string> _discrete_parameters_RB;
  std::vector<Real> _discrete_parameter_values_in_RB;
  std::map<std::string,std::vector<Real>> _discrete_parameter_values_RB;

  Real _rel_training_tolerance_RB;
  Real _abs_training_tolerance_RB;
  std::vector<Real> _continuous_parameter_min_values_RB;
  std::vector<Real> _continuous_parameter_max_values_RB;

  unsigned int _qp;
  unsigned int _i;
  unsigned int _j;
  unsigned int _num_elems;

  MooseVariable & _var;
  const VariableTestValue & _test;

  const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system;
};

#endif //DWARFELEPHANTCOMPUTEEIMINNERPRODUCTMATRIXSTEADYSTATE_H
