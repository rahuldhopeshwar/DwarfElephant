#include "RBNodalBC.h"
#include "MooseVariable.h"
#include "Assembly.h"

template<>
InputParameters validParams<RBNodalBC>()
{
  InputParameters params = validParams<NodalBC>();

  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");
  params.addRequiredParam<FunctionName>("cache_stiffness_matrix", "");

  return params;
}


RBNodalBC::RBNodalBC(const InputParameters & parameters) :
    NodalBC(parameters),
    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject")),
    _function(&getFunction("cache_stiffness_matrix"))
{

    _cache_stiffness_matrix = dynamic_cast<CacheStiffnessMatrix *>(_function);
}

void
RBNodalBC::computeResidual(NumericVector<Number> & residual)
{
  if (_var.isNodalDefined())
  {
    dof_id_type & dof_idx = _var.nodalDofIndex();
    _qp = 0;
    Real res = computeQpResidual();
    residual.set(dof_idx, res);


    if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
    {
      _cache_stiffness_matrix->cacheResidual(dof_idx, res);
    }

    if (_has_save_in)
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i=0; i<_save_in.size(); i++)
        _save_in[i]->sys().solution().set(_save_in[i]->nodalDofIndex(), res);
    }
  }
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

    if(_initialize_rb_system._offline_stage)
    {
      if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
      {
       const std::vector<BoundaryName> & _boundary_names = boundaryNames();
       _cache_stiffness_matrix->resizeSubdomainCaches(_initialize_rb_system._qa);

        for(unsigned int _i = 0; _i != _boundary_names.size(); _i++)
        {
          if (_boundary_names[_i] == "bottom")
          {
            _cache_stiffness_matrix->cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val, 0);
          }
          else if (_boundary_names[_i] == "top")
          {
            _cache_stiffness_matrix->cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val, _initialize_rb_system._qa-1);
          }
//
//          else
//          {
//            const std::set<SubdomainID> & _node_block_ids = _mesh.getNodeBlockIds(*_current_node);
//
//	        for(std::set<SubdomainID>::const_iterator it = _node_block_ids.begin();
//                it != _node_block_ids.end(); it++)
//            {
//                 _initialize_rb_system._jacobian_subdomain[*it] ->close();
//                 _initialize_rb_system._jacobian_subdomain[*it] -> zero_rows(cached_rows, cached_val);
//            }
//          }
       }
//        _cache_stiffness_matrix -> cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val);
      }
    }


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
