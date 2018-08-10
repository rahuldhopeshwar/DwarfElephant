#include "DwarfElephantLiftingFunction.h"

template <>
InputParameters
validParams<DwarfElephantLiftingFunction>()
{
  InputParameters params = validParams<Function>();
  params.addClassDescription("This class implements the lifitng function for variable z-values.");
  params.addRequiredParam<FunctionName>("boundary_function", "The function that describes the value distribution of the boundary.");
  params.addRequiredParam<unsigned int>("ID_reference_layer", "The ID of the reference layer.");
  params.addParam<unsigned int>("ID_data_layer", 0,"The ID of the data layer.");
  params.addParam<std::vector<unsigned int>>("dim_data_array", "The dimensions of the data set.");
  params.addParam<std::vector<Real>>("step_sizes", "The step size of the data set in x and y direction.");
  params.addParam<Real>("tolerance", 0.1,"Tolerance for the search procedure.");
  params.addParam<bool>("access_multiple_times", false, "Whether the data needs to be accessed multiple times.");
  return params;
}

DwarfElephantLiftingFunction::DwarfElephantLiftingFunction(const InputParameters & parameters)
  : Function(parameters),
    FunctionInterface(this),
    _ID_reference_layer(getParam<unsigned int>("ID_reference_layer")),
    _ID_data_layer(getParam<unsigned int>("ID_data_layer")),
    _data_array_dimensions(getParam<std::vector<unsigned int>>("dim_data_array")),
    _step_sizes(getParam<std::vector<Real>>("step_sizes")),
    _tolerance(getParam<Real>("tolerance")),
    _access_multiple_times(getParam<bool>("access_multiple_times"))
{
  const FunctionName & _boundary_function_name = getParam<FunctionName>("boundary_function");
  _boundary_function = &getFunctionByName(_boundary_function_name);

  const std::vector<dof_id_type> * _reference_node_ids = &_sc_fe_problem.mesh().getNodeList(_ID_reference_layer);

  for(unsigned int i=0;  i< _reference_node_ids[0].size(); i++)
  {
    Node & _reference_node = _sc_fe_problem.mesh().nodeRef(_reference_node_ids[0][i]);
    _x_coord_reference_layer.push_back(_reference_node(0));
    _y_coord_reference_layer.push_back(_reference_node(1));
    _z_coord_reference_layer.push_back(_reference_node(2));
  }

  if(_access_multiple_times)
  {
    _data_array.resize(_data_array_dimensions[0]);
    for(unsigned int i=0; i<_data_array_dimensions[0]; i++)
      _data_array[i].resize(_data_array_dimensions[1]);

    const std::vector<dof_id_type> * _reference_node_ids = &_sc_fe_problem.mesh().getNodeList(_ID_data_layer);

    for(unsigned int i=0;  i< _reference_node_ids[0].size(); i++)
    {
      Node & _reference_node = _sc_fe_problem.mesh().nodeRef(_reference_node_ids[0][i]);
      Real _value = interpolateZRef(_reference_node(0), _reference_node(1));
      Real _x = std::round(_reference_node(0)/_step_sizes[0]);
      Real _y = std::round(_reference_node(1)/_step_sizes[1]);
      _data_array[_x][_y]=_value;
    }
  }

}

void
DwarfElephantLiftingFunction::initialSetup()
{
}

Real
DwarfElephantLiftingFunction::value(Real /*t*/, const Point & p)
{
  Real _x = std::round(p(0)/_step_sizes[0]);
  Real _y = std::round(p(1)/_step_sizes[1]);
  Real _z_ref = _data_array[_x][_y];

  // for(unsigned int i=0; i<_x_coord_reference_layer.size(); i++)
  //   if ((std::fabs(_x_coord_reference_layer[i] - p(0)) < _tolerance) && (std::fabs(_y_coord_reference_layer[i] - p(1)) < _tolerance))
  //     _z_ref = _z_coord_reference_layer[i];

  Real _value = _boundary_function->value(0,p)*((p(2)-_z_ref)/std::fabs(p(2)-_z_ref));
  return _value;
}

RealGradient
DwarfElephantLiftingFunction::gradient(Real /*t*/, const Point & p)
{
  Real _x = std::round(p(0)/_step_sizes[0]);
  Real _y = std::round(p(1)/_step_sizes[1]);
  Real _z_ref = _data_array[_x][_y];

  Real _z_value = (1/std::fabs(_z_ref-p(2)))-((pow((p(2)-_z_ref),2))/(pow((std::fabs(_z_ref-p(2))),3)));
  RealGradient _f_gradient(0,0,_z_value);
  RealGradient _gradient = /*(_boundary_function->gradient(0,p)*((p(2)-_z_ref)/std::fabs(p(2)-_z_ref))) +*/
               _boundary_function->value(0,p) * _f_gradient;
  return _gradient;
}

Real
DwarfElephantLiftingFunction::interpolateZRef(Real _x_coord, Real _y_coord)
{
  Real sum = 0.0;
  Real values = 0.0;
  Real distance = 0.0;
  std::vector<Real> weights;

  for(unsigned int i=0; i<_x_coord_reference_layer.size(); i++)
  {
    if ((std::fabs(_x_coord_reference_layer[i] - _x_coord) < _tolerance) && (std::fabs(_y_coord_reference_layer[i] - _y_coord) < _tolerance))
      return _z_coord_reference_layer[i];

    distance = pow((_x_coord_reference_layer[i] - _x_coord),2)+pow((_y_coord_reference_layer[i] - _y_coord),2);
    weights.push_back(1.0 / pow(distance,2));
  }

  for (unsigned int i = 0; i < weights.size(); ++i)
  {
    sum += weights[i];
    values += weights[i] * _z_coord_reference_layer[i];
  }
  return values / sum;
}
