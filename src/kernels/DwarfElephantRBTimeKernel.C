#include "DwarfElephantRBTimeKernel.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariable.h"
#include "SystemBase.h"

// libMesh includes
#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<DwarfElephantRBTimeKernel>()
{
  InputParameters params = validParams<DwarfElephantRBKernel>();
  return params;
}

DwarfElephantRBTimeKernel::DwarfElephantRBTimeKernel(const InputParameters & parameters) : DwarfElephantRBKernel(parameters) {}

void
DwarfElephantRBTimeKernel::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number(), Moose::KT_TIME);
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
    for (unsigned int i = 0; i < _save_in.size(); i++)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}
