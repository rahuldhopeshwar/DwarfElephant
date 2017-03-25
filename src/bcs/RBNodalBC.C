#include "RBNodalBC.h"
#include "MooseVariable.h"
#include "Assembly.h"

template<>
InputParameters validParams<RBNodalBC>()
{
  InputParameters params = validParams<NodalBC>();
  
  params += validParams<BlockRestrictable>();
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");

  return params;
}


RBNodalBC::RBNodalBC(const InputParameters & parameters) :
    NodalBC(parameters),
    BlockRestrictable(parameters),
    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject"))
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
    
    if(_initialize_rb_system._offline_stage)
    {
      if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
      {
        const std::set<SubdomainID> & _node_block_ids = _mesh.getNodeBlockIds(*_current_node);
	
	for(std::set<SubdomainID>::const_iterator it = _node_block_ids.begin();
            it != _node_block_ids.end(); it++)
	{
	  //for(unsigned _q = 0; _q < _initialize_rb_system._qa; _q++)
	    //if (*it == _q)
              _initialize_rb_system._jacobian_subdomain[*it] -> set(cached_row, cached_row, cached_val);
        }
	
	_initialize_rb_system._inner_product_matrix -> set(cached_row, cached_row, cached_val);
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
