#include "DwarfElephantLiftingFunction.h"
#include "DwarfElephantBoundaryConditionFileReader.h"

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
  params.addParam<std::string>("dx_distance_file", "Name and path of the dx file of the distance.");
  params.addParam<std::string>("dy_distance_file", "Name and path of the dy file of the distance.");
  params.addParam<std::string>("dx_file", "Name and path of the dx file.");
  params.addParam<std::string>("dy_file", "Name and path of the dy file.");
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
    _dx_distance_file(getParam<std::string>("dx_distance_file")),
    _dy_distance_file(getParam<std::string>("dy_distance_file")),
    _dx_file(getParam<std::string>("dx_file")),
    _dy_file(getParam<std::string>("dy_file")),
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

  _dx = fileParserGradients(_dx_file);
  _dy = fileParserGradients(_dy_file);

  _dx_distance = fileParserGradients(_dx_distance_file);
  _dy_distance = fileParserGradients(_dy_distance_file);

  if(_access_multiple_times)
  {
    _data_array.resize(_data_array_dimensions[0]);
    _distance_data_array.resize(_data_array_dimensions[0]);
    _dx_distance_data_array.resize(_data_array_dimensions[0]);
    _dy_distance_data_array.resize(_data_array_dimensions[0]);
    _dx_data_array.resize(_data_array_dimensions[0]);
    _dy_data_array.resize(_data_array_dimensions[0]);

    for(unsigned int i=0; i<_data_array_dimensions[0]; i++)
    {
      _data_array[i].resize(_data_array_dimensions[1]);
      _distance_data_array[i].resize(_data_array_dimensions[1]);
      _dx_distance_data_array[i].resize(_data_array_dimensions[1]);
      _dy_distance_data_array[i].resize(_data_array_dimensions[1]);
      _dx_data_array[i].resize(_data_array_dimensions[1]);
      _dy_data_array[i].resize(_data_array_dimensions[1]);
    }

    const std::vector<dof_id_type> * _reference_node_ids = &_sc_fe_problem.mesh().getNodeList(_ID_data_layer);

    for(unsigned int i=0;  i< _reference_node_ids[0].size(); i++)
    {
      Node & _reference_node = _sc_fe_problem.mesh().nodeRef(_reference_node_ids[0][i]);
      Real z_ref = interpolateZRef(_reference_node(0), _reference_node(1));
      _data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=z_ref;
      _distance_data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=std::fabs(_reference_node(2)-z_ref);
      _dx_distance_data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=findGradientDistance(_reference_node(0), _reference_node(1))[0];
      _dy_distance_data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=findGradientDistance(_reference_node(0), _reference_node(1))[1];
      _dx_data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=findGradientDepth(_reference_node(0), _reference_node(1))[0];
      _dy_data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=findGradientDepth(_reference_node(0), _reference_node(1))[1];
    }
  }
}

void
DwarfElephantLiftingFunction::initialSetup()
{
  // std::ofstream myfile("depth_data_array.dat");
  // myfile << "i " << "j " << "distance" << std::endl;
  //
  // for(unsigned int j = 0; j< _data_array_dimensions[1]; j++)
  //   for(unsigned int i = 0; i< _data_array_dimensions[0]; i++)
  //     myfile << i << " " << j << " " <<_data_array[i][j] << std::endl;
}

Real
DwarfElephantLiftingFunction::value(Real /*t*/, const Point & p)
{
  Real _x = std::round(p(0)/_step_sizes[0]);
  Real _y = std::round(p(1)/_step_sizes[1]);
  Real _z_ref = _data_array[_x][_y];
  Real _distance = _distance_data_array[_x][_y];

  Real _value = _boundary_function->value(0,p)*((p(2)-_z_ref)/_distance);
  return _value;
}

RealGradient
DwarfElephantLiftingFunction::gradient(Real /*t*/, const Point & p)
{
  Real _x = std::round(p(0)/_step_sizes[0]);
  Real _y = std::round(p(1)/_step_sizes[1]);
  Real _z_ref = _data_array[_x][_y];
  Real _distance = _distance_data_array[_x][_y];
  Real dx_distance = _dx_distance_data_array[_x][_y];
  Real dy_distance = _dy_distance_data_array[_x][_y];
  Real dx = _dx_data_array[_x][_y];
  Real dy = _dy_data_array[_x][_y];
  Real _fx = (dx*_distance)+((p(2)-_z_ref)*dx_distance);
  Real _fy = (dy*_distance)+((p(2)-_z_ref)*dy_distance);

  RealGradient _f_gradient(_fx, _fy, (-1.0 * _distance));
  RealGradient _gradient = _boundary_function->gradient(0,p)*((p(2)-_z_ref)/_distance) -
               _boundary_function->value(0,p) * 1.0/pow(_distance,2) * _f_gradient;
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

std::vector<Real>
DwarfElephantLiftingFunction::findGradientDistance(Real _x_coord, Real _y_coord)
{
  std::vector<Real> _gradient;
  _gradient.resize(2);
  for (unsigned int i=0; i < _x_coord_reference_layer.size(); ++i)
  {
    if ((std::fabs(_x_coord_reference_layer[i] - _x_coord) < _tolerance) && (std::fabs(_y_coord_reference_layer[i] - _y_coord) < _tolerance))
    {
      _gradient[0]=_dx_distance[i];
      _gradient[1]=_dy_distance[i];
      return _gradient;
    }
  }
   mooseError("Gradient for Point(", _x_coord, ", ", _y_coord, ") ",
              "could not be matched.");
  return _gradient;
}

std::vector<Real>
DwarfElephantLiftingFunction::findGradientDepth(Real _x_coord, Real _y_coord)
{
  std::vector<Real> _gradient;
  _gradient.resize(2);
  for (unsigned int i=0; i < _x_coord_reference_layer.size(); ++i)
  {
    if ((std::fabs(_x_coord_reference_layer[i] - _x_coord) < _tolerance) && (std::fabs(_y_coord_reference_layer[i] - _y_coord) < _tolerance))
    {
      _gradient[0]=_dx[i];
      _gradient[1]=_dy[i];
      return _gradient;
    }
  }
   mooseError("Gradient for Point(", _x_coord, ", ", _y_coord, ") ",
              "could not be matched.");
  return _gradient;
}

std::vector<Real>
DwarfElephantLiftingFunction::fileParserGradients(std::string & file)
{
  // define all necessary parameters
  std::string _line;
  std::vector<Real> _gradient_data;

  // Check inputFile and create input stream.
  MooseUtils::checkFileReadable(file);
  std::ifstream _input_file;

  // Open the inputFile and check it.
  _input_file.open(file.c_str());
  if (!_input_file.good())
    mooseError("Error while opening the file '", file, "' within the DwarfElephantFileReader function.");

  // Read the file.
  while (std::getline(_input_file, _line))
  {
    if (_line[0] != '#')
    {
      // Parse the data to the corresponding vectors
      std::istringstream _iss(_line);
      Real _data_value;
      while (_iss >> _data_value)
        _gradient_data.push_back(_data_value);
    }
  }
  _input_file.close();

  return _gradient_data;
}
