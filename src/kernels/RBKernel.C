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

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
RBKernel::RBKernel(const InputParameters & parameters) :
    Kernel(parameters),
    _use_displaced(getParam<bool>("use_displaced")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
    _block_ids(this->blockIDs()),
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

  // Add the calcualted matrices to the Aq matrices from the RB system.
  _initialize_rb_system._jacobian_subdomain[*_block_ids.begin()] -> add_matrix(_local_ke, _var.dofIndices());
  _initialize_rb_system._jacobian_subdomain[*_block_ids.begin()] -> close();

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
