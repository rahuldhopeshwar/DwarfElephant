/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBKernel inherits from the Kernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * stiffness matrix is needed separated in its subdomain contributions. In
 * addition it overwrites the function computeResidual.
 */

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/threads.h"
#include "libmesh/quadrature.h"

////MOOSE includes
#include "Assembly.h"
#include "MooseVariable.h"
#include "Problem.h"
#include "SubProblem.h"
#include "SystemBase.h"

//MOOSE includes (DwarfElephant package)
#include "RBKernel.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<RBKernel>()
{
  InputParameters params = validParams<Kernel>();

  params.addClassDescription("Overwrites the function computeJacobian. This is required because for the RB method the stiffness matrix needs to be saved in its subdomain contributions.");
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");
  params.addParam<unsigned int>("ID_Aq", 0, "ID of the current stiffness matrix");
  params.addParam<unsigned int>("ID_Fq", 0, "ID of the current stiffness matrix");
  params.addParam<bool>("matrix_separation_according_to_subdomains", true, "Tells whether the stiffness matrix is separated according to the subdomain_ids");
  params.addParam<bool>("vector_separation_according_to_subdomains", true, "Tells whether the load vector is separated according to the subdomain_ids");
  params.addRequiredParam<Real>("max_x","Maximum extension of the volume of interest in x-direction.");
  params.addRequiredParam<Real>("min_x","Minimum extension of the volume of interest in x-direction.");
  params.addRequiredParam<Real>("max_y","Maximum extension of the volume of interest in y-direction.");
  params.addRequiredParam<Real>("min_y","Minimum extension of the volume of interest in y-direction.");
  params.addRequiredParam<Real>("max_z","Maximum extension of the volume of interest in z-direction.");
  params.addRequiredParam<Real>("min_z","Minimum extension of the volume of interest in z-direction.");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
RBKernel::RBKernel(const InputParameters & parameters) :
    Kernel(parameters),
    _use_displaced(getParam<bool>("use_displaced")),
    _matrix_separation_according_to_subdomains(getParam<bool>("matrix_separation_according_to_subdomains")),
    _vector_separation_according_to_subdomains(getParam<bool>("vector_separation_according_to_subdomains")),
    _ID_Aq(getParam<unsigned int>("ID_Aq")),
    _ID_Fq(getParam<unsigned int>("ID_Fq")),
    _max_x(getParam<Real>("max_x")),
    _min_x(getParam<Real>("min_x")),
    _max_y(getParam<Real>("max_y")),
    _min_y(getParam<Real>("min_y")),
    _max_z(getParam<Real>("max_z")),
    _min_z(getParam<Real>("min_z")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
    _block_ids(this->blockIDs()),
    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject"))

{
}

///-------------------------------------------------------------------------
void
RBKernel::initialSetup()
{
  // Error messages
  if (_block_ids.size()>1)
  {
      mooseError("For the RB method the stiffness matrix has to be saved separatly for each subdomain. Therefore each RBKernel and each inheriting Kernel needs to be defined individually for each block. You defined the Kernel for more than one block, please change your specifications in the Input file.");
  }

  if(_initialize_rb_system._exec_flags[0] != EXEC_INITIAL)
    mooseError("The initialization of the RB system has to be executed on 'initial'. "
               "You defined a wrong state in your 'execute_on' line in the Inputfile. "
               "Please, correct your settings.");

  // Defining the IDs of the stiffness matrix and load vectors in case of subdomain separation
  if(_matrix_separation_according_to_subdomains)
    _ID_Aq = *_block_ids.begin();

  if(_vector_separation_according_to_subdomains)
    _ID_Fq = *_block_ids.begin();

  _output_volume = (_max_x - _min_x) * (_max_y - _min_y) * (_max_z - _min_z);
}

void
RBKernel::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  _local_out.resize(re.size());
  _local_out.zero();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); _i++)
  {
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    {
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();
      _local_out(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual() / _output_volume;
    }
  }

  re += _local_re;

  if(_initialize_rb_system._offline_stage)
  {
//     Real _output_volume = (_max_x - _min_x) * (_max_y - _min_y) * (_max_z - _min_z);
//    _fe_problem.mesh();
    // Add the calculated matrices to the Aq matrices from the RB system.
    if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
    {
        _initialize_rb_system._residuals[_ID_Fq] -> add_vector(_local_re, _var.dofIndices());
        _initialize_rb_system._outputs[0] -> add_vector(_local_out, _var.dofIndices());
    }
  }

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (const auto & var : _save_in)
      var->sys().solution().add_vector(_local_re, var->dofIndices());
  }
}

void
RBKernel::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  precalculateJacobian();
  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _phi.size(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();

  ke += _local_ke;

  if(_initialize_rb_system._offline_stage)
  {
    // Add the calculated matrices to the Aq matrices from the RB system.
    if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
    {
        _initialize_rb_system._jacobian_subdomain[_ID_Aq] -> add_matrix(_local_ke, _var.dofIndices());
    }
  }

 if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i=0; i<rows; i++)
      diag(i) = _local_ke(i,i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (const auto & var : _diag_save_in)
      var->sys().solution().add_vector(diag, var->dofIndices());
  }
}
///----------------------------------PDEs-----------------------------------
// For the PDEs zero is implemented, since this Kernel shall be used for any
// RB problem. The problem specific PDEs are implemented in separate Kernels.

Real
RBKernel::computeQpJacobian()
{
  return 0;
}

Real
RBKernel::computeQpResidual()
{
  return 0;
}
