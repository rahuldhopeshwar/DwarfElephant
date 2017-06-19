#include "RBIntegratedBC.h"
#include "SubProblem.h"
#include "SystemBase.h"
#include "MooseVariable.h"
#include "Assembly.h"

// libMesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<RBIntegratedBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addParam<unsigned int>("ID_Aq", 0, "ID of the current stiffness matrix");
  params.addParam<unsigned int>("ID_Fq", 0, "ID of the current stiffness matrix");

  return params;
}

RBIntegratedBC::RBIntegratedBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _use_displaced(getParam<bool>("use_displaced")),
    _simulation_type(getParam<std::string>("simulation_type")),
    _ID_Aq(getParam<unsigned int>("ID_Aq")),
    _ID_Fq(getParam<unsigned int>("ID_Fq")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es())
{
}

RBIntegratedBC::~RBIntegratedBC()
{
}

void
RBIntegratedBC::initialSetup()
{
  _output_volume = 0.00777;
}

void
RBIntegratedBC::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  _local_out.resize(re.size());
  _local_out.zero();

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    for (_i = 0; _i < _test.size(); _i++)
      _local_re(_i) += _JxW[_qp]*_coord[_qp]*computeQpResidual();

  re += _local_re;

//  if ((_min_x <= _centroid(0)) && (_centroid(0) <= _max_x) &&
//      (_min_y <= _centroid(1)) && (_centroid(1) <= _max_y))
    for (_i = 0; _i < _test.size(); _i++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _local_out(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual() / _output_volume;


  if(_simulation_type == "steady")  // SteadyState
  {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
      // Add the calculated vectors to the vectors from the RB system.
      if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        _initialize_rb_system._residuals[_ID_Fq] -> add_vector(_local_re, _var.dofIndices());
//        _initialize_rb_system._residuals[0] -> add_vector(_local_re, _var.dofIndices());
//        _initialize_rb_system._outputs[0] -> add_vector(_local_out, _var.dofIndices());
  }

  else if (_simulation_type == "transient") // Transient
  {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
      // Add the calculated vectors to the vectors from the RB system.
      if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        _initialize_rb_system._residuals[_ID_Fq] -> add_vector(_local_re, _var.dofIndices());
//        _initialize_rb_system._residuals[0] -> add_vector(_local_re, _var.dofIndices());
//        _initialize_rb_system._outputs[0] -> add_vector(_local_out, _var.dofIndices());
  }


  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i=0; i<_save_in.size(); i++)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}

void
RBIntegratedBC::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        _local_ke(_i, _j) += _JxW[_qp]*_coord[_qp]*computeQpJacobian();

  ke += _local_ke;

  if(_simulation_type == "steady")  // Steady State
  {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
    // Add the calculated matrices to the Aq matrices from the RB system.
    if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
        _initialize_rb_system._jacobian_subdomain[_ID_Aq] -> add_matrix(_local_ke, _var.dofIndices());
   }

  else if(_simulation_type == "transient") // Transient
  {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
    // Add the calculated matrices to the Aq matrices from the RB system.
    if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
        _initialize_rb_system._jacobian_subdomain[_ID_Aq] -> add_matrix(_local_ke, _var.dofIndices());
  }

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i=0; i<rows; i++)
      diag(i) = _local_ke(i,i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i=0; i<_diag_save_in.size(); i++)
      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }
}

void
RBIntegratedBC::computeJacobianBlock(unsigned int jvar)
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
    for (_i=0; _i<_test.size(); _i++)
      for (_j=0; _j<_phi.size(); _j++)
      {
        if (_var.number() == jvar)
          ke(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpJacobian();
        else
          ke(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpOffDiagJacobian(jvar);
      }
}

void
RBIntegratedBC::computeJacobianBlockScalar(unsigned int jvar)
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);

  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < jv.order(); _j++)
        ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpOffDiagJacobian(jvar);
}

Real
RBIntegratedBC::computeQpResidual()
{
  return 0;
}

Real
RBIntegratedBC::computeQpJacobian()
{
  return 0;
}

Real
RBIntegratedBC::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
  return 0;
}
