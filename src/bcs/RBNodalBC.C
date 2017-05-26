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
  params.addRequiredParam<FunctionName>("cache_boundaries", "");
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("mesh_modified", false, "Determine whether you operate on a modified mesh or not.");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
RBNodalBC::RBNodalBC(const InputParameters & parameters) :
    NodalBC(parameters),
    _mesh_modified(getParam<bool>("mesh_modified")),
    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject")),
    _function(&getFunction("cache_boundaries"))
{

    _cache_boundaries = dynamic_cast<CacheBoundaries *>(_function);
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

    if(_initialize_rb_system._offline_stage)
    {
      if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
      {

       _cache_boundaries->resizeSubdomainVectorCaches(_initialize_rb_system._qf);
       _cache_boundaries->cacheResidual(dof_idx, -res);

      if(_mesh_modified)
      {
        // MOOSE modified mesh
        const std::vector<BoundaryName> & _boundary_names = boundaryNames();

        for(unsigned int _i = 0; _i != _boundary_names.size(); _i++)
        {
          if (_boundary_names[_i] == "bottom")
            _cache_boundaries->cacheSubdomainResidual(dof_idx, -res, 0);

          else if (_boundary_names[_i] == "top")
            _cache_boundaries->cacheSubdomainResidual(dof_idx, -res, _initialize_rb_system._qf-1);
        }
       }

       else
       {
         // external Mesh
         const std::set< SubdomainID > & _node_boundary_list = _mesh.getNodeBlockIds(*_current_node);
         _cache_boundaries->cacheSubdomainResidual(dof_idx, -res, *_node_boundary_list.begin());
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

    if(_initialize_rb_system._offline_stage)
    {
      if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
      {
       _cache_boundaries -> cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val);

       _cache_boundaries->resizeSubdomainMatrixCaches(_initialize_rb_system._qa);

        if (_mesh_modified)
        {
          // MOOSE modified mesh
          const std::vector<BoundaryName> & _boundary_names = boundaryNames();

          for(unsigned int _i = 0; _i != _boundary_names.size(); _i++)
          {
            if (_boundary_names[_i] == "bottom")
            {
              _cache_boundaries->cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val, 0);
            }
            else if (_boundary_names[_i] == "top")
            {
              _cache_boundaries->cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val, _initialize_rb_system._qa-1);
            }
          }
        }

        else
        {
          // external mesh
          const std::set< SubdomainID > & _node_boundary_list = _mesh.getNodeBlockIds(*_current_node);
          _cache_boundaries->cacheSubdomainStiffnessMatrixContribution(cached_row, cached_row, cached_val,*_node_boundary_list.begin());
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
