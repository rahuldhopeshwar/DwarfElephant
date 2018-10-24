#include "DwarfElephantElementalVariableValuesAlongLine.h"

// MOOSE includes
#include "LineSegment.h"
#include "RayTracing.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

registerMooseObject("DwarfElephantApp", DwarfElephantElementalVariableValuesAlongLine);

template <>
InputParameters
validParams<DwarfElephantElementalVariableValuesAlongLine>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addRequiredParam<VariableName>("variable", "The variable you want to be monitor.");
  params.addRequiredParam<Point>("start", "The beginning of the line");
  params.addRequiredParam<Point>("end", "The end of the line");

  return params;
}

DwarfElephantElementalVariableValuesAlongLine::DwarfElephantElementalVariableValuesAlongLine(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _elem_ids(declareVector("elem_ids")),
    _elem_x(declareVector("elem_x")),
    _elem_y(declareVector("elem_y")),
    _elem_z(declareVector("elem_z")),
    _elem_values(declareVector("elem_values")),
    _var_name(parameters.get<VariableName>("variable")),
    _start(getParam<Point>("start")),
    _end(getParam<Point>("end"))
{
}

void
DwarfElephantElementalVariableValuesAlongLine::initialize()
{
  _elem_ids.clear();
  _elem_x.clear();
  _elem_y.clear();
  _elem_z.clear();
  _elem_values.clear();
}

void
DwarfElephantElementalVariableValuesAlongLine::execute()
{
  Real value = 0;
  std::vector<Elem *> intersected_elems;
  std::vector<LineSegment> segments;

  std::unique_ptr<PointLocatorBase> pl = _fe_problem.mesh().getPointLocator();
  Moose::elementsIntersectedByLine(
      _start, _end, _fe_problem.mesh(), *pl, intersected_elems, segments);

  unsigned int num_elems = intersected_elems.size();

  _elem_ids.resize(num_elems);
  _elem_x.resize(num_elems);
  _elem_y.resize(num_elems);
  _elem_z.resize(num_elems);
  _elem_values.resize(num_elems);

  for (unsigned int i = 0; i < num_elems; i++)
  {
    // Get Element Values
    _subproblem.prepare(intersected_elems[i], _tid);
    _subproblem.reinitElem(intersected_elems[i], _tid);

    // In case your are using a MOOSE version older than Feb 20, 2018
    // use the following line
    //MooseVariable & var = _subproblem.getVariable(_tid, _var_name);
    MooseVariable & var = _subproblem.getStandardVariable(_tid, _var_name);
    const VariableValue & u = var.sln();

    unsigned int n = u.size();
    for (unsigned int i = 0; i < n; i++)
      value += u[i];
    value /= n;

    gatherSum(value);

    _elem_ids[i] = intersected_elems[i]->id();
    _elem_x[i] = intersected_elems[i]->centroid().operator()(0);
    _elem_y[i] = intersected_elems[i]->centroid().operator()(1);
    _elem_z[i] = intersected_elems[i]->centroid().operator()(2);
    _elem_values[i] = value;
  }
}
