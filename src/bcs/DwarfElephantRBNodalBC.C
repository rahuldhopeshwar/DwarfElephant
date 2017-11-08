/**
 * This BC is required to use the RB method as it is provided by the
 * RB libMesh package. The RBNodalBC inherits from the NodalBC class. It
 * overwrites the function computeJacobian because for the RB method the
 * stiffness matrix is needed separated in its subdomain contributions. In
 * addition it overwrites the function computeResidual.
 */

///---------------------------------INCLUDES--------------------------------
#include "DwarfElephantRBNodalBC.h"
#include "MooseVariable.h"
#include "Assembly.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBNodalBC>()
{
  InputParameters params = validParams<NodalBC>();

  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addParam<unsigned int>("ID_Fq", 0 , "ID if the load vector.");
  params.addParam<unsigned int>("ID_Aq", 0 , "ID if the stiffness matrix.");
  params.addParam<unsigned int>("ID_Mq", 0 , "ID if the mass matrix.");
  params.addParam<unsigned int>("ID_Oq", 0, "ID of the current output vector");
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("matrix_seperation_according_to_subdomains", true, "Tells whether the stiffness matrix is separated according to the subdomain_ids");
  params.addParam<bool>("compute_output",false,"Determines whether an output function is used or not");
  params.addParam<Real>("max_x", 0.,"Maximum extension of the volume of interest in x-direction.");
  params.addParam<Real>("min_x", 0.,"Minimum extension of the volume of interest in x-direction.");
  params.addParam<Real>("max_y", 0.,"Maximum extension of the volume of interest in y-direction.");
  params.addParam<Real>("min_y", 0.,"Minimum extension of the volume of interest in y-direction.");
  params.addParam<Real>("max_z", 0.,"Maximum extension of the volume of interest in z-direction.");
  params.addParam<Real>("min_z", 0.,"Minimum extension of the volume of interest in z-direction.");


  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantRBNodalBC::DwarfElephantRBNodalBC(const InputParameters & parameters) :
    NodalBC(parameters),
    _matrix_seperation_according_to_subdomains(getParam<bool>("matrix_seperation_according_to_subdomains")),
    _compute_output(getParam<bool>("compute_output")),
    _simulation_type(getParam<std::string>("simulation_type")),
    _ID_Fq(getParam<unsigned int>("ID_Fq")),
    _ID_Aq(getParam<unsigned int>("ID_Aq")),
    _ID_Mq(getParam<unsigned int>("ID_Mq")),
    _ID_Oq(getParam<unsigned int>("ID_Oq")),
    _max_x(getParam<Real>("max_x")),
    _min_x(getParam<Real>("min_x")),
    _max_y(getParam<Real>("max_y")),
    _min_y(getParam<Real>("min_y")),
    _max_z(getParam<Real>("max_z")),
    _min_z(getParam<Real>("min_z"))
{
    _rb_problem = cast_ptr<DwarfElephantRBProblem *>(&_fe_problem);
}

///-------------------------------------------------------------------------
void
DwarfElephantRBNodalBC::computeResidual(NumericVector<Number> & residual)
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

//          _rb_problem->rbAssembly(_ID_Fq).cacheResidual(dof_idx, -res);
          _rb_problem->rbAssembly().cacheResidual(dof_idx, -res, _ID_Fq);
        }
      }
    }

    else if (_simulation_type == "transient")
    {
      const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");

      if(_initialize_rb_system._offline_stage)
        if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        {}
//          _rb_problem->rbAssembly(_ID_Fq).cacheResidual(dof_idx, -res);
    }

    if (_has_save_in)
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i=0; i<_save_in.size(); i++)
        _save_in[i]->sys().solution().set(_save_in[i]->nodalDofIndex(), res);
    }
  }

  if(_compute_output)
    computeOutput();
}

void
DwarfElephantRBNodalBC::computeOutput()
{

  if ((_min_x <= _current_node->operator()(0)) && (_current_node->operator()(0) <= _max_x) &&
      (_min_y <= _current_node->operator()(1)) && (_current_node->operator()(2) <= _max_y) &&
      (_min_z <= _current_node->operator()(2)) && (_current_node->operator()(2) <= _max_z))
  {
    if (_var.isNodalDefined())
    {
      dof_id_type & dof_idx = _var.nodalDofIndex();
      _qp = 0;
      Real res = computeQpResidual();

      if (_simulation_type == "steady")  // Steady State
      {
        const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");

        if(_initialize_rb_system._offline_stage)
          if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
          {}
//          _rb_problem->rbAssembly(_ID_Oq).cacheOutput(_var.nodalDofIndex(), -res);
	  //_rb_problem->rbAssembly(_ID_Oq).cacheOutput(dof_idx, -res);
      }

      else if (_simulation_type == "transient")
      {
        const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");

        if(_initialize_rb_system._offline_stage)
          if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
          {}
//            _rb_problem->rbAssembly(_ID_Oq).cacheOutput(dof_idx, -res);
      }
    }
  }
}

void
DwarfElephantRBNodalBC::computeJacobian()
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

          if (_matrix_seperation_according_to_subdomains)
          {
          const std::set< SubdomainID > & _node_boundary_list = _mesh.getNodeBlockIds(*_current_node);
          for (std::set<SubdomainID>::const_iterator it = _node_boundary_list.begin();
             it != _node_boundary_list.end(); ++it)
            _rb_problem->rbAssembly().cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val, *it);
//            _rb_problem->rbAssembly(*it).cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val);
           }
           else
	        _rb_problem->rbAssembly().cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val, _ID_Aq);
//	        _rb_problem->rbAssembly(_ID_Aq).cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val);
        }
      }
    }

    else if (_simulation_type == "transient")
    {
      const DwarfElephantInitializeRBSystemTransient& _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");

      if(_initialize_rb_system._offline_stage)
        if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
        {
//          _rb_problem->rbAssembly(_ID_Aq).cacheStiffnessMatrixContribution(cached_row, cached_row, cached_val);
//          _rb_problem->rbAssembly(_ID_Mq).cacheMassMatrixContribution(cached_row, cached_row, cached_val);
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
DwarfElephantRBNodalBC::computeQpJacobian()
{
  return 1.;
}

Real
DwarfElephantRBNodalBC::computeQpResidual()
{
  return 0.;
}
