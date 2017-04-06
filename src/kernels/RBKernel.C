/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBKernel inherits from the Kernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * stiffness matrix is needed separated in its subdomain contributions.
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
  params.addRequiredParam<unsigned int>("subdomain", "The active subdomain");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
RBKernel::RBKernel(const InputParameters & parameters) :
    Kernel(parameters),
    _use_displaced(getParam<bool>("use_displaced")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
    _block_ids(this->blockIDs()),
    _block(getParam<unsigned int>("subdomain")),
    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject"))

{
}

///-------------------------------------------------------------------------
void
RBKernel::initialSetup()
{
  if (_block_ids.size()>1)
  {
      mooseError("For the RB method the stiffness matrix has to be saved separatly for each subdomain. Therefore each RBKernel and each inheriting Kernel needs to be defined individually for each block. You defined the Kernel for more than one block, please change your specifications in the Input file.");
  }

  if(_initialize_rb_system._exec_flags[0] != EXEC_INITIAL)
    mooseError("The initialization of the RB system has to be executed on 'initial'. "
               "You defined a wrong state in your 'execute_on' line in the Inputfile. "
               "Please, correct your settings.");
}

void
RBKernel::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); _i++)
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();

  re += _local_re;

  if(_initialize_rb_system._offline_stage)
  {
    // Add the calculated matrices to the Aq matrices from the RB system.
    if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
    {
////        _initialize_rb_system._residuals[*_block_ids.begin()];
        _initialize_rb_system._residuals[*_block_ids.begin()] -> add_vector(_local_re, _var.dofIndices());
    }
      if (_initialize_rb_system._compliant)
      {
////        _initialize_rb_system._outputs[*_block_ids.begin()];
        _initialize_rb_system._outputs[*_block_ids.begin()] -> add_vector(_local_re, _var.dofIndices());
      }
//
//      else if (!_initialize_rb_system._compliant)
//        mooseError ("Currently, the implementation handles only one output term and only the compliant case.");
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
        _initialize_rb_system._jacobian_subdomain[*_block_ids.begin()] -> add_matrix(_local_ke, _var.dofIndices());
        _initialize_rb_system._inner_product_matrix -> add_matrix(_local_ke, _var.dofIndices());
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
