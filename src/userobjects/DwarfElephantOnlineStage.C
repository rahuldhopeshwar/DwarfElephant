 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantOnlineStage.h"

template<>
InputParameters validParams<DwarfElephantOnlineStage>()
{
  InputParameters params = validParams<NodalUserObject>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");
  params.addRequiredParam<NonlinearVariableName>("variable", "");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");

  return params;
}

DwarfElephantOnlineStage::DwarfElephantOnlineStage(const InputParameters & params):
  NodalUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _system_name(getParam<std::string>("system")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _mesh_ptr(&_fe_problem.mesh()),
  _nodal_bcs(false),
  _variable_name(params.get<NonlinearVariableName>("variable")),
  _var(_fe_problem.getVariable(_tid, _variable_name)),
  _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject"))
{
}

void
DwarfElephantOnlineStage::onlineStage()
{
  ////    _rb_eval.legacy_read_offline_data_from_files();
////    RBParameters _online_mu;
////
////    _online_mu.set_value("mu0", _online_mu0_parameters);
////    _rb_eval.set_parameters(_online_mu);
////    _rb_eval.rb_solve(_online_N);
////
////    _rb_eval.print_parameters();
}

void
DwarfElephantOnlineStage::initialize()
{
}

void
DwarfElephantOnlineStage::execute()
{
  unsigned int _off_diag_jacobian_value = 0;
  _console << _current_node->id() << std::endl;

  NodeRange & _nodes = *_mesh_ptr->getActiveNodeRange();

  for (const auto & _node : _nodes)
    if (_current_node->id() != _node->id())
    {
      _initialize_rb_system._inner_product_matrix->set(_current_node->id(), _node->id(), _off_diag_jacobian_value);
    }
}

void
DwarfElephantOnlineStage::threadJoin(const UserObject & y)
{
}

void
DwarfElephantOnlineStage::finalize()
{
}
