#include "DwarfElephantRBPointValue.h"

// MOOSE includes
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SubProblem.h"

#include "libmesh/system.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBPointValue);

template <>
InputParameters
validParams<DwarfElephantRBPointValue>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<Point>("point", "The physical point where we monitor");
  params.addRequiredParam<VariableName>(
      "variable", "The name of the variable that this postprocessor operates on.");
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
    _var_number(_subproblem.getVariable(_tid, parameters.get<VariableName>("variable")).number()),
    _outputid(getParam<unsigned int>("outputid")),
    _simulation_type(getParam<std::string>("simulation_type"))
{
}

void
DwarfElephantRBPointValue::execute()
{
  if(_simulation_type=="steady") {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    assignPoint(_initialize_rb_system.getOutputs(), _outputid, _point);
  } else {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    assignPoint(_initialize_rb_system.getOutputs(), _outputid, _point);
  }
}

void
DwarfElephantRBPointValue::assignPoint(const std::vector<std::vector<NumericVector <Number> *> > _outputs,
  unsigned int _outputid, Point _point, bool _insistOnSuccess)
{
  if (_outputid >= _outputs.size())
    mooseError("The outputid you defined is out of range.");

  // This function must be called on every processor because it is not known in
  // advance which processors holds the point of interest.
  parallel_object_only();

  //It is important to check that every processor agrees on the point
  #ifndef NDEBUG
  libmesh_assert(comm().verify(_point(0)));
  #if LIBMESH_DIM > 1
  libmesh_assert(comm().verify(_point(1)));
  #endif
  #if LIBMESH_DIM > 2
  libmesh_assert(comm().verify(_point(2)));
  #endif
  #endif // NDEBUG

  // Use an existing PointLocator or create a new one.
  std::unique_ptr<PointLocatorBase> _locator_ptr = _mesh.getMesh().sub_point_locator();
  PointLocatorBase & _locator = *_locator_ptr;

  if (!_insistOnSuccess || !_mesh.getMesh().is_serial())
    _locator.enable_out_of_mesh_mode();

  // Get a pointer to the element that contains the point.
  const Elem * _element = _locator(_point);

  // Get the DofMap for the element.
  const DofMap & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();
  std::vector<dof_id_type> dof_indices;

  dof_map.dof_indices(_element, dof_indices, _var_number);

  // Prepare the calculation of the weights.
  FEType fe_type = dof_map.variable_type(_var_number);

  Point coor = FEInterface::inverse_map(_element->dim(), fe_type, _element, _point);
  FEComputeData fe_data(_fe_problem.es(), coor);
  FEInterface::compute_data(_element->dim(), fe_type, _element, fe_data);

  // Assign the correct weights to the output of interest vector.
  if (_element && (_element->processor_id() == processor_id()))
  {
    for(unsigned int _q=0; _q < _outputs[_outputid].size(); _q++)
    {
      for(unsigned int _n=0; _n < dof_indices.size(); _n++)
      {
        _outputs[_outputid][_q]->set(dof_indices[_n],fe_data.shape[_n]);
      }
      // _outputs[_outputid][_q]->close();
    }
  }
}
