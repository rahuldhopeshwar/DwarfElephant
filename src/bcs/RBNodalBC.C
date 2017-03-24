#include "RBNodalBC.h"
#include "MooseVariable.h"
#include "Assembly.h"

template<>
InputParameters validParams<RBNodalBC>()
{
  InputParameters params = validParams<NodalBC>();

  return params;
}


RBNodalBC::RBNodalBC(const InputParameters & parameters) :
    NodalBC(parameters)
{
}

void
RBNodalBC::computeJacobian()
{
  // We call the user's computeQpJacobian() function and store the
  // results in the _assembly object. We can't store them directly in
  // the element stiffness matrix, as they will only be inserted after
  // all the assembly is done.
  if (_var.isNodalDefined())
  {
    _qp = 0;
    Real cached_val = computeQpJacobian();
    dof_id_type cached_row = _var.nodalDofIndex();

    // Cache the user's computeQpJacobian() value for later use.
    _fe_problem.assembly(0).cacheJacobianContribution(cached_row, cached_row, cached_val);

    if (_has_diag_save_in)
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i=0; i<_diag_save_in.size(); i++)
        _diag_save_in[i]->sys().solution().set(_diag_save_in[i]->nodalDofIndex(), cached_val);
    }
  }
}


Real
RBNodalBC::computeQpJacobian()
{
  return 1.;
}

Real
RBNodalBC::computeQpResidual()
{
  return 0.;
}
