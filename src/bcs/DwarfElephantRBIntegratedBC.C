#include "DwarfElephantRBIntegratedBC.h"
#include "SubProblem.h"
#include "SystemBase.h"
#include "MooseVariable.h"
#include "Assembly.h"

// libMesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<DwarfElephantRBIntegratedBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params += validParams<BlockRestrictable>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addParam<unsigned int>("ID_Aq", 0, "ID of the current stiffness matrix");
  params.addParam<unsigned int>("ID_Mq", 0, "ID of the current mass matrix");
  params.addParam<unsigned int>("ID_Fq", 0, "ID of the current load vector");
  params.addParam<unsigned int>("ID_Oq", 0, "ID of the current output vector");
  params.addParam<bool>("matrix_seperation_according_to_subdomains", true, "Tells whether the stiffness matrix is separated according to the subdomain_ids");
  params.addParam<bool>("compute_output",false,"Determines whether an output function is used or not");
  params.addParam<Real>("max_x", 0.,"Maximum extension of the volume of interest in x-direction.");
  params.addParam<Real>("min_x", 0.,"Minimum extension of the volume of interest in x-direction.");
  params.addParam<Real>("max_y", 0.,"Maximum extension of the volume of interest in y-direction.");
  params.addParam<Real>("min_y", 0.,"Minimum extension of the volume of interest in y-direction.");
  params.addParam<Real>("max_z", 0.,"Maximum extension of the volume of interest in z-direction.");
  params.addParam<Real>("min_z", 0.,"Minimum extension of the volume of interest in z-direction.");

  return params;
}

DwarfElephantRBIntegratedBC::DwarfElephantRBIntegratedBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    BlockRestrictable(parameters),
    _use_displaced(getParam<bool>("use_displaced")),
    _matrix_seperation_according_to_subdomains(getParam<bool>("matrix_seperation_according_to_subdomains")),
    _compute_output(getParam<bool>("compute_output")),
    _simulation_type(getParam<std::string>("simulation_type")),
    _ID_first_block(*_fe_problem.mesh().meshSubdomains().begin()),
    _ID_Aq(getParam<unsigned int>("ID_Aq")),
    _ID_Mq(getParam<unsigned int>("ID_Mq")),
    _ID_Fq(getParam<unsigned int>("ID_Fq")),
    _ID_Oq(getParam<unsigned int>("ID_Oq")),
    _max_x(getParam<Real>("max_x")),
    _min_x(getParam<Real>("min_x")),
    _max_y(getParam<Real>("max_y")),
    _min_y(getParam<Real>("min_y")),
    _max_z(getParam<Real>("max_z")),
    _min_z(getParam<Real>("min_z")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es())
{
}

DwarfElephantRBIntegratedBC::~DwarfElephantRBIntegratedBC()
{
}

void
DwarfElephantRBIntegratedBC::initialSetup()
{
  if(_compute_output)
    _output_volume = (_max_x - _min_x) * (_max_y - _min_y) * (_max_z - _min_z);
}

void
DwarfElephantRBIntegratedBC::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    for (_i = 0; _i < _test.size(); _i++)
      _local_re(_i) += _JxW[_qp]*_coord[_qp]*computeQpResidual();

  re += _local_re;

  if(_simulation_type == "steady")  // SteadyState
  {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
      // Add the calculated vectors to the vectors from the RB system.
      if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        _initialize_rb_system._residuals[_ID_Fq] -> add_vector(_local_re, _var.dofIndices());

  }

  else if (_simulation_type == "transient") // Transient
  {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
      // Add the calculated vectors to the vectors from the RB system.
      if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        _initialize_rb_system._residuals[_ID_Fq] -> add_vector(_local_re, _var.dofIndices());
  }


  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i=0; i<_save_in.size(); i++)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }

 // if(_compute_output)
 //   computeOutput();
}

void
DwarfElephantRBIntegratedBC::computeOutput()
{
  DenseVector<Number> & out = _assembly.residualBlock(_var.number());
  _local_out.resize(out.size());
  _local_out.zero();

  Point _centroid = _current_elem->centroid();

  if ((_min_x <= _centroid(0)) && (_centroid(0) <= _max_x) &&
      (_min_y <= _centroid(1)) && (_centroid(1) <= _max_y) &&
      (_min_z <= _centroid(2)) && (_centroid(2) <= _max_z))
  {
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      for (_i = 0; _i < _test.size(); _i++)
        _local_out(_i) += _JxW[_qp]*_coord[_qp]*computeQpResidual();

    out += _local_out;
  }

  if(_simulation_type == "steady")  // SteadyState
  {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
      // Add the calculated vectors to the vectors from the RB system.
      if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
      {
        _initialize_rb_system._outputs[0][_ID_Oq] -> add_vector(_local_out, _var.dofIndices());
//        _console << _local_out << std::endl;
     }
  }

  else if (_simulation_type == "transient") // Transient
  {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
      // Add the calculated vectors to the vectors from the RB system.
      if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
         _initialize_rb_system._outputs[0][_ID_Oq] -> add_vector(_local_out, _var.dofIndices());
  }
}

void
DwarfElephantRBIntegratedBC::computeJacobian()
{
  if(_matrix_seperation_according_to_subdomains)
    _ID_Aq = _current_elem->subdomain_id() - _ID_first_block;

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
    {
        _initialize_rb_system._jacobian_subdomain[_ID_Aq] -> add_matrix(_local_ke, _var.dofIndices());
        _initialize_rb_system._mass_matrix_subdomain[_ID_Mq] -> add_matrix(_local_ke, _var.dofIndices());
    }
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
DwarfElephantRBIntegratedBC::computeJacobianBlock(unsigned int jvar)
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
DwarfElephantRBIntegratedBC::computeJacobianBlockScalar(unsigned int jvar)
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);

  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < jv.order(); _j++)
        ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpOffDiagJacobian(jvar);
}

Real
DwarfElephantRBIntegratedBC::computeQpResidual()
{
  return 0;
}

Real
DwarfElephantRBIntegratedBC::computeQpJacobian()
{
  return 0;
}

Real
DwarfElephantRBIntegratedBC::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
  return 0;
}
