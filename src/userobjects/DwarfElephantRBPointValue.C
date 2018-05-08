#include "DwarfElephantRBPointValue.h"

// MOOSE includes
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SubProblem.h"

#include "libmesh/system.h"

template <>
InputParameters
validParams<DwarfElephantRBPointValue>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<Point>("point", "The physical point where we monitor");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");
  params.addParam<unsigned int>("outputid",0,"The ID of the output of interest.");
  params.addClassDescription("Outputs an elemental variable value at a particular location");

  return params;
}

DwarfElephantRBPointValue::DwarfElephantRBPointValue(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _mesh(_subproblem.mesh()),
    _point(getParam<Point>("point")),
    _outputid(getParam<unsigned int>("outputid")),
    _simulation_type(getParam<std::string>("simulation_type"))
{
}

void
DwarfElephantRBPointValue::execute()
{
  if(_simulation_type=="steady") {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    assignPoint(_initialize_rb_system.getOutputs());
  } else {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    assignPoint(_initialize_rb_system.getOutputs());
  }
}

void
DwarfElephantRBPointValue::assignPoint(const std::vector<std::vector<NumericVector <Number> *> > _outputs)
{
  // if (_elementid.size() != _outputs.size())
  //   mooseError("The number of elements to extract does not match the number of outputs used in the redution.");
  //
  // for(unsigned int i=0; i < _outputs.size(); i++)
  // {
  //   Elem * _element = _mesh.getMesh().query_elem_ptr(_elementid[i]);
  //   const Node *const* node = _element->get_nodes();
  //
  //   if (_element && (_element->processor_id() == processor_id()))
  //   {
  //
  //     Real _value = 1./_element->n_nodes();
  //
  //     for(unsigned int _q=0; _q < _outputs[i].size(); _q++)
  //     {
  //       for(unsigned int _n=0; _n < _element->n_nodes(); _n++)
  //       {
  //         _outputs[i][_q]->set(node[_n]->id(),_value);
  //         _outputs[i][_q]->close();
  //       }
  //     }
  //   }
  // }
}
