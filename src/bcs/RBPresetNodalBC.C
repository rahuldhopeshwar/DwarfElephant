#include "RBPresetNodalBC.h"

#include "libmesh/numeric_vector.h"

template<>
InputParameters validParams<RBPresetNodalBC>()
{
  InputParameters p = validParams<RBNodalBC>();
  return p;
}


RBPresetNodalBC::RBPresetNodalBC(const InputParameters & parameters) :
  RBNodalBC(parameters)
{

}

void
RBPresetNodalBC::computeValue(NumericVector<Number> & current_solution)
{
  dof_id_type & dof_idx = _var.nodalDofIndex();
  _qp = 0;
  current_solution.set(dof_idx, computeQpValue());
}

Real
RBPresetNodalBC::computeQpResidual()
{
  return _u[_qp] - computeQpValue();
}

Real
RBPresetNodalBC::computeQpValue()
{
  return 0;
}
