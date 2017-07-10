/**
 * This BC is required to use the RB method as it is provided by the
 * RB libMesh package. The RBNodalBC inherits from the NodalBC class. It
 * overwrites the function computeJacobian because for the RB method the
 * stiffness matrix is needed separated in its subdomain contributions. In
 * addition it overwrites the function computeResidual.
 */

///---------------------------------INCLUDES--------------------------------
#include "RBNodalBC.h"
#include "MooseVariable.h"
#include "Assembly.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<RBNodalBC>()
{
  InputParameters params = validParams<NodalBC>();

  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");
//  params.addRequiredParam<FunctionName>("cache_boundaries", "");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addParam<unsigned int>("ID_Fq", 0 , "ID if the load vector.");
  params.addParam<unsigned int>("ID_Aq", 0 , "ID if the stiffness matrix.");
  params.addParam<unsigned int>("ID_Mq", 0 , "ID if the mass matrix.");
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
RBNodalBC::RBNodalBC(const InputParameters & parameters) :
    NodalBC(parameters),
    _simulation_type(getParam<std::string>("simulation_type")),
    _ID_Fq(getParam<unsigned int>("ID_Fq")),
    _ID_Aq(getParam<unsigned int>("ID_Aq")),
    _ID_Mq(getParam<unsigned int>("ID_Mq"))
//    _function(&getFunction("cache_boundaries"))
{

//    _cache_boundaries = dynamic_cast<CacheBoundaries *>(_function);
    _rb_problem = cast_ptr<DwarfElephantRBProblem *>(&_fe_problem);

    //_rb_problem->newRBAssemblyArray(_fe_problem.getNonlinearSystemBase());

}

///-------------------------------------------------------------------------
void
RBNodalBC::computeResidual(NumericVector<Number> & residual)
{
  if (_var.isNodalDefined())
  {
    dof_id_type & dof_idx = _var.nodalDofIndex();
    _qp = 0;
    Real res = computeQpResidual();
    residual.set(dof_idx, res);

    if (_simulation_type == "steady")  // Steady State
    {
      const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");

      if(_initialize_rb_system._offline_stage)
      {
        if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        {

          _rb_problem->rbAssembly(_ID_Fq).cacheResidual(dof_idx, -res);
	  //_initialize_rb_system._residuals[_ID_Fq]->set(dof_idx,-res);

         // external Mesh
//         const std::set< SubdomainID > & _node_boundary_list = _mesh.getNodeBlockIds(*_current_node);
//         for (std::set<SubdomainID>::const_iterator it = _node_boundary_list.begin();
//              it != _node_boundary_list.end(); ++it)
//           _cache_boundaries->cacheSubdomainResidual(dof_idx, -res, *it - 0);
//           _cache_boundaries->cacheSubdomainResidual(dof_idx, -res, *it - _ID_first_block);
//           _cache_boundaries->cacheSubdomainResidual(dof_idx, -res, _ID_Fq);
        }
      }
    }

    else if (_simulation_type == "transient")
    {
      const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");

      if(_initialize_rb_system._offline_stage)
      {
        if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        {
          //_cache_boundaries->resizeSubdomainVectorCaches(_initialize_rb_system._qf);
//          _cache_boundaries->cacheResidual(dof_idx, -res);
//          _cache_boundaries->cacheSubdomainResidual(dof_idx, -res, _ID_Fq);
        }
      }
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

    if (_simulation_type == "steady")
    {
      const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");

      if(_initialize_rb_system._offline_stage)
      {
        if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0 )
        {

	      _rb_problem->rbAssembly(_ID_Aq).cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val);

          // external mesh
//          const std::set< SubdomainID > & _node_boundary_list = _mesh.getNodeBlockIds(*_current_node);
//          for (std::set<SubdomainID>::const_iterator it = _node_boundary_list.begin();
//               it != _node_boundary_list.end(); ++it)
//            _cache_boundaries->cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val,*it - _ID_first_block);
//            _cache_boundaries->cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val, _ID_Aq);    _rb_problem->rbAssembly(0).cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val, _ID_Aq);
        }
      }
    }

    else if (_simulation_type == "transient")
    {
      const DwarfElephantInitializeRBSystemTransient& _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");

      if(_initialize_rb_system._offline_stage)
      {
        if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
        {
//          _cache_boundaries -> cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val);
          //_cache_boundaries->resizeSubdomainStiffnessMatrixCaches(_initialize_rb_system._qa);
          //_cache_boundaries->resizeSubdomainMassMatrixCaches(_initialize_rb_system._qm);
          //_cache_boundaries->cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val, _ID_Aq);
//          _cache_boundaries->cacheSubdomainMassMatrixContribution(cached_row, cached_row, cached_val, _ID_Mq);
        }
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
