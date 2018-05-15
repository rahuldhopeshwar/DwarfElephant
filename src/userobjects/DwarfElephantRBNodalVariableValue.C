#include "DwarfElephantRBNodalVariableValue.h"

// MOOSE includes
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SubProblem.h"

#include "libmesh/node.h"

template <>
InputParameters
validParams<DwarfElephantRBNodalVariableValue>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<std::vector<unsigned int>>("nodeid", "The ID of the node where we monitor");
  params.addParam<Real>("scale_factor", 1, "A scale factor to be applied to the variable");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");
  params.addClassDescription("Outputs values of a nodal variable at a particular location");

  return params;
}

DwarfElephantRBNodalVariableValue::DwarfElephantRBNodalVariableValue(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _mesh(_subproblem.mesh()),
    _nodeid(getParam<std::vector<unsigned int>>("nodeid")),
    _scale_factor(getParam<Real>("scale_factor")),
    _simulation_type(getParam<std::string>("simulation_type"))
{
  // This class may be too dangerous to use if renumbering is enabled,
  // as the nodeid parameter obviously depends on a particular
  // numbering.
    if (_mesh.getMesh().allow_renumbering())
      mooseError("DwarfElephantRBNodalVariableValue should only be used when node renumbering is disabled.");

  for (unsigned int i = 0; i < _nodeid.size(); i++)
  {
    Node * _node_ptr = _mesh.getMesh().query_node_ptr(_nodeid[i]);

    bool found_node_ptr = _node_ptr;
    _communicator.max(found_node_ptr);

    if (!found_node_ptr)
      mooseError("Node #",
                 getParam<unsigned int>("nodeid"),
                 " specified in '",
                 name(),
                 "' not found in the mesh!");
  }
}

void
DwarfElephantRBNodalVariableValue::execute()
{
  if(_simulation_type=="steady") {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    assignNode(_initialize_rb_system.getOutputs());
  } else {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    assignNode(_initialize_rb_system.getOutputs());
  }
}

void
DwarfElephantRBNodalVariableValue::assignNode(const std::vector<std::vector<NumericVector <Number> *> > _outputs)
{
  if (_nodeid.size() != _outputs.size())
    mooseError("The number of nodes to extract does not match the number of outputs used in the redution.");

  Real _value = 1 * _scale_factor;

  for(unsigned int i=0; i < _outputs.size(); i++)
  {
    Node * _node_ptr = _mesh.getMesh().query_node_ptr(_nodeid[i]);

    if (_node_ptr && _node_ptr->processor_id() == processor_id())
    {
      for(unsigned int _q=0; _q < _outputs[i].size(); _q++)
      {
        _outputs[i][_q]->set(_nodeid[i],_value);
        // _outputs[i][_q]->close();
      }
    }
  }
}
