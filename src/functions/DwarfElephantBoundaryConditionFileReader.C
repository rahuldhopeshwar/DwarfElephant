#include "DwarfElephantBoundaryConditionFileReader.h"

template <>
InputParameters
validParams<DwarfElephantBoundaryConditionFileReader>()
{
  InputParameters params = validParams<Function>();
  params.addClassDescription("This file reader is required for reading in the boundary conditions.");
  params.addRequiredParam<std::string>("file", "Name and path of the data file, valid delimiters are space and tab-space.");
  params.addParam<int>("dimension", 2, "Spatial dimension of the data set." );
  // params.addParam<std::string>("delimiter_to_replace", "In case you data file contains a not tab separated data set.");
  return params;
}

DwarfElephantBoundaryConditionFileReader::DwarfElephantBoundaryConditionFileReader(const InputParameters & parameters)
  : Function(parameters),
  _file(getParam<std::string>("file")),
  _dimension(getParam<int>("dimension"))
  // _delimiter_to_replace(getParam<std::string>("delimiter_to_replace"))
{
  if(_dimension<1 && _dimension > 3)
    mooseError("The spatial dimension are incorrect.");
  fileParser();
}

Real
DwarfElephantBoundaryConditionFileReader::value(Real /*t*/, const Point & p)
{
  return findValue(p(0), p(1)/*, p(2)*/);
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

      if(_dimension==1)
        _x_coordinates.push_back(_all_data[0]);
      else if (_dimension == 2)
      {
        _x_coordinates.push_back(_all_data[0]);
        _y_coordinates.push_back(_all_data[1]);
      }
      // else
      // {
      //   _x_coordinates.push_back(_all_data[0]);
      //   _y_coordinates.push_back(_all_data[1]);
      //   _z_coordinates.push_back(_all_data[2]);
      // }
      _variable_values.push_back(_all_data[_dimension]);
    }
  }
  _input_file.close();
  _num_points = _x_coordinates.size();
}

Real
DwarfElephantBoundaryConditionFileReader::findValue(Real _x_coord, Real _y_coord /*, Real _z_coord*/)
{
  for (unsigned int i=0; i < _num_points; ++i)
  {
    if ((std::fabs(_x_coordinates[i] - _x_coord) < 0.1) && (std::fabs(_y_coordinates[i] - _y_coord) < 0.1) /*&& (std::fabs(_z_coordinates[i] - _z_coord) < 0.1)*/)
      return _variable_values[i];
  }
   mooseError("Point could not be matched.");
  return 0;
}
