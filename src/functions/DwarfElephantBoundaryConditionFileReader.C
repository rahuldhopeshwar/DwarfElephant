#include "DwarfElephantBoundaryConditionFileReader.h"

registerMooseObject("DwarfElephantApp", DwarfElephantBoundaryConditionFileReader);

template <>
InputParameters
validParams<DwarfElephantBoundaryConditionFileReader>()
{
  InputParameters params = validParams<Function>();
  params.addClassDescription("This file reader is required for reading in the boundary conditions.");
  params.addRequiredParam<std::string>("file", "Name and path of the data file, valid delimiters are space and tab-space.");
  params.addParam<std::string>("dx_file", "","Name and path of the dx file.");
  params.addParam<std::string>("dy_file", "","Name and path of the dx file.");
  params.addParam<unsigned int>("dimension", 2, "Spatial dimension of the data set." );
  params.addParam<unsigned int>("ID_data_layer", 0,"The ID of the data layer.");
  params.addParam<std::vector<unsigned int>>("dim_data_array", "The dimensions of the data set.");
  params.addParam<std::vector<Real>>("step_sizes", "The step size of the data set in x and y direction.");
  params.addParam<Real>("tolerance", 0.1,"Tolerance for the search procedure.");
  params.addParam<bool>("interpolate", false, "If true than the data will be interpolated in space.");
  params.addParam<bool>("access_multiple_times", false, "Whether the data needs to be accessed multiple times.");
  params.addParam<bool>("gradients", false, "Whether the gradients have to be read.");
  return params;
}

DwarfElephantBoundaryConditionFileReader::DwarfElephantBoundaryConditionFileReader(const InputParameters & parameters)
  : Function(parameters),
  _file(getParam<std::string>("file")),
  _dx_file(getParam<std::string>("dx_file")),
  _dy_file(getParam<std::string>("dy_file")),
  _dimension(getParam<unsigned int>("dimension")),
  _ID_data_layer(getParam<unsigned int>("ID_data_layer")),
  _data_array_dimensions(getParam<std::vector<unsigned int>>("dim_data_array")),
  _step_sizes(getParam<std::vector<Real>>("step_sizes")),
  _tolerance(getParam<Real>("tolerance")),
  _interpolate(getParam<bool>("interpolate")),
  _access_multiple_times(getParam<bool>("access_multiple_times")),
  _gradients(getParam<bool>("gradients"))
{
  if(_dimension!=2)
    mooseError("The spatial dimension are incorrect.");

  fileParser();

  if(_gradients)
  {
    _dx = fileParserGradients(_dx_file);
    _dy = fileParserGradients(_dy_file);
  }

    if(_gradients && !_access_multiple_times)
      mooseError("The reading process of the gradients is currently only supported",
                 " for the multiple access option.");

  if(_access_multiple_times)
  {
    if(_data_array_dimensions.size()!=2)
      mooseError("Wrong dimensions of the data array.");

    if(_step_sizes.size()!=2)
      mooseError("Wrong dimensions of the step sizes of the data array.");

    _data_array.resize(_data_array_dimensions[0]);

    if(_gradients)
    {
      _dx_data_array.resize(_data_array_dimensions[0]);
      _dy_data_array.resize(_data_array_dimensions[0]);
    }

    for(unsigned int i=0; i<_data_array_dimensions[0]; i++)
    {
      _data_array[i].resize(_data_array_dimensions[1]);

      if(_gradients)
      {
        _dx_data_array[i].resize(_data_array_dimensions[1]);
        _dy_data_array[i].resize(_data_array_dimensions[1]);
      }
    }

    const std::vector<dof_id_type> * _reference_node_ids = &_sc_fe_problem.mesh().getNodeList(_ID_data_layer);


    for(unsigned int i=0;  i< _reference_node_ids[0].size(); i++)
    {
      Node & _reference_node = _sc_fe_problem.mesh().nodeRef(_reference_node_ids[0][i]);
      _data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=findValue(_reference_node(0), _reference_node(1));

      if(_gradients)
      {
        _dx_data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=findGradient(_reference_node(0), _reference_node(1))[0];
        _dy_data_array[std::round(_reference_node(0)/_step_sizes[0])][std::round(_reference_node(1)/_step_sizes[1])]=findGradient(_reference_node(0), _reference_node(1))[1];
      }
    }
  }
}

Real
DwarfElephantBoundaryConditionFileReader::value(Real /*t*/, const Point & p) const
{
  if(!_interpolate && !_access_multiple_times)
    return findValue(p(0), p(1));
  else if(_interpolate && !_access_multiple_times)
    return interpolateValue(p(0), p(1));
  else
    return _data_array[std::round(p(0)/_step_sizes[0])][std::round(p(1)/_step_sizes[1])];
}

RealGradient
DwarfElephantBoundaryConditionFileReader::gradient(Real /*t*/, const Point & p) const
{
  Real _x = std::round(p(0)/_step_sizes[0]);
  Real _y = std::round(p(1)/_step_sizes[1]);
  if(_access_multiple_times)
  {
     return (_dx_data_array[_x][_y], _dy_data_array[_x][_y], 0);
  }
  return (0,0,0);
}

void
DwarfElephantBoundaryConditionFileReader::fileParser()
{
  // define all necessary parameters
  std::string _line;

  // Check inputFile and create input stream.
  MooseUtils::checkFileReadable(_file);
  std::ifstream _input_file;

  // Open the inputFile and check it.
  _input_file.open(_file.c_str());
  if (!_input_file.good())
    mooseError("Error while opening the file '", _file, "' within the DwarfElephantFileReader function.");

  // Read the file.
  while (std::getline(_input_file, _line))
  {
    if (_line[0] != '#')
    {
      // Parse the data to the corresponding vectors
      std::istringstream _iss(_line);
      std::vector<Real> _all_data;
      Real _data_value;
      while (_iss >> _data_value)
        _all_data.push_back(_data_value);

      _x_coordinates.push_back(_all_data[0]);
      _y_coordinates.push_back(_all_data[1]);
      _variable_values.push_back(_all_data[3]);
    }
  }
  _input_file.close();
  _num_points = _x_coordinates.size();
}

std::vector<Real>
DwarfElephantBoundaryConditionFileReader::fileParserGradients(std::string & file)
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
  _num_points = _x_coordinates.size();

  return _gradient_data;
}

Real
DwarfElephantBoundaryConditionFileReader::findValue(Real _x_coord, Real _y_coord) const
{
  for (unsigned int i=0; i < _num_points; ++i)
  {
    if ((std::fabs(_x_coordinates[i] - _x_coord) < _tolerance) && (std::fabs(_y_coordinates[i] - _y_coord) < _tolerance))
      return _variable_values[i];
  }
   mooseError("Point(", _x_coord, ", ", _y_coord, ") ",
              "could not be matched.");
  return 0;
}

std::vector<Real>
DwarfElephantBoundaryConditionFileReader::findGradient(Real _x_coord, Real _y_coord) const
{
  std::vector<Real> _gradient;
  _gradient.resize(2);
  for (unsigned int i=0; i < _num_points; ++i)
  {
    if ((std::fabs(_x_coordinates[i] - _x_coord) < _tolerance) && (std::fabs(_y_coordinates[i] - _y_coord) < _tolerance))
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

Real
DwarfElephantBoundaryConditionFileReader::interpolateValue(Real _x_coord, Real _y_coord) const
{
  Real sum = 0.0;
  Real values = 0.0;
  Real distance = 0.0;
  std::vector<Real> weights;

  for (unsigned int i=0; i < _num_points; ++i)
  {
    if ((std::fabs(_x_coordinates[i] - _x_coord) < _tolerance) && (std::fabs(_y_coordinates[i] - _y_coord) < _tolerance))
      return _variable_values[i];

    distance = pow((_x_coordinates[i] - _x_coord),2)+pow((_y_coordinates[i] - _y_coord),2);
    weights.push_back(1.0 / pow(distance,2));
  }

  for (unsigned int i = 0; i < weights.size(); ++i)
  {
    sum += weights[i];
    values += weights[i] * _variable_values[i];
  }
  return values / sum;
}
