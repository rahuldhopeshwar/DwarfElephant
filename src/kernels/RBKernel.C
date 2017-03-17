/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBKernel inherits from the Kernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * whole stiffness matrix is needed and not only the diagonal entries.
 */

///---------------------------------INCLUDES--------------------------------
//MOOSE includes (DwarfElephant package)
#include "RBKernel.h"

////MOOSE includes
#include "Assembly.h"
#include "MooseVariable.h"
#include "Problem.h"
#include "SubProblem.h"
#include "SystemBase.h"

// libMesh includes
#include "libmesh/threads.h"
#include "libmesh/quadrature.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<RBKernel>()
{
  InputParameters params = validParams<Kernel>();

  params.addClassDescription("Overwrites the function computeJacobian. This is required because for the RB method the entire stiffness matrix needs to be saved and not only the diagonal entries.");
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
RBKernel::RBKernel(const InputParameters & parameters) :
    Kernel(parameters),
    _system_name(getParam<std::string>("system")),
    _use_displaced(getParam<bool>("use_displaced")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
    _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name))

{
}

///-------------------------------------------------------------------------
void
RBKernel::timestepSetup()
{
   _rb_con = &_es.get_system<DwarfElephantRBConstruction>("RBSystem");
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

  if (_has_diag_save_in)
  {
//    unsigned int rows = ke.m();
//    unsigned int columns = ke.n();
//    DenseVector<Number> actualRow(rows);
//
//    for (unsigned int j=0; j<columns; j++)
//    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (const auto & var : _diag_save_in)
      {
        _jacobian = _rb_con->get_Aq(0);
        _jacobian -> add_matrix(_local_ke, var->dofIndices());
        _jacobian->close();
        _console << *_jacobian << std::endl;
//        _rb_con->get_Aq(0)->add_matrix(_local_ke, var->dofIndices());
//        for (unsigned int i=0; i<rows; i++)
//        {
//          actualRow(i) = _local_ke(i,j);
//        }
//       var->sys().solution().add_vector(actualRow, var->dofIndices());

       }
//     }
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
