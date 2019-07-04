#include "DwarfElephantSingleCoordinateFileReader.h"

registerMooseObject("DwarfElephantApp", DwarfElephantSingleCoordinateFileReader);

template <>
InputParameters
validParams<DwarfElephantSingleCoordinateFileReader>()
{
  InputParameters params = validParams<Function>();
  params.addClassDescription("This file reader is required for reading in the boundary conditions.");
  params.addRequiredParam<std::string>("file", "Name and path of the data file, valid delimiters are space and tab-space.");
  params.addParam<unsigned int>("coordinate", 2, "The coordinate that should be read from file.");
  params.addParam<Real>("tolerance", 0.1,"Tolerance for the search procedure.");
  params.addParam<bool>("interpolate", false, "If true than the data will be interpolated in space.");
  params.addParam<bool>("access_multiple_times", false, "Whether the data needs to be accessed multiple times.");
  params.addParam<bool>("gradients", false, "Whether the gradients have to be read.");
  return params;
}

DwarfElephantSingleCoordinateFileReader::DwarfElephantSingleCoordinateFileReader(const InputParameters & parameters)
  : Function(parameters),
  _file(getParam<std::string>("file")),
  _coordinate(getParam<unsigned int>("coordinate")),
  _tolerance(getParam<Real>("tolerance")),
  _interpolate(getParam<bool>("interpolate"))
{
  fileParser();
}

Real
DwarfElephantSingleCoordinateFileReader::value(Real /*t*/, const Point & p) const
{
  if(!_interpolate)
    return findValue(p(_coordinate));
  else
    return interpolateValue(p(_coordinate));
}

void
DwarfElephantSingleCoordinateFileReader::fileParser()
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
      // std::vector<Real> _all_data;
      Real _data_value;
      while (_iss >> _data_value)
        _coordinates.push_back(_data_value);
        // _all_data.push_back(_data_value);

      // _coordinates.push_back(_all_data[0]);
    }
  }
  _input_file.close();
  _num_points = _coordinates.size();
}

Real
DwarfElephantSingleCoordinateFileReader::findValue(Real _coord) const
{
  for (unsigned int i=0; i < _num_points; ++i)
  {
    if ((std::fabs(_coordinates[i] - _coord) < _tolerance))
      return _coordinates[i];
  }
   mooseError("Point(", _coord, ") ",
              "could not be matched.");
  return 0;
}

Real
DwarfElephantSingleCoordinateFileReader::interpolateValue(Real _coord) const
{
  Real sum = 0.0;
  Real values = 0.0;
  Real distance = 0.0;
  std::vector<Real> weights;

  for (unsigned int i=0; i < _num_points; ++i)
  {
    if ((std::fabs(_coordinates[i] - _coord) < _tolerance))
      return _coordinates[i];

    distance = pow((_coordinates[i] - _coord),2);
    weights.push_back(1.0 / pow(distance,2));
  }

  for (unsigned int i = 0; i < weights.size(); ++i)
  {
    sum += weights[i];
    values += weights[i] * _coordinates[i];
  }
  return values / sum;
}
