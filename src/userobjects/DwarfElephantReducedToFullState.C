 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantReducedToFullState.h"

template<>
InputParameters validParams<DwarfElephantReducedToFullState>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<std::string>("file_base", "The reduced state vector file without extension.");
  params.addParam<bool>("read_binary_data",true,"Determines if the format of the state file is binary or ascii.");
  return params;
}

DwarfElephantReducedToFullState::DwarfElephantReducedToFullState(const InputParameters & params):
  GeneralUserObject(params),
  _file_base(getParam<std::string>("file_base")),
  _read_binary_data(getParam<bool>("read_binary_data"))
{
}

void
DwarfElephantReducedToFullState::initialize()
{
  read_reduced_state(_file_base, _read_binary_data);
}

void
DwarfElephantReducedToFullState::execute()
{
}

void
DwarfElephantReducedToFullState::finalize()
{
}

void
DwarfElephantReducedToFullState::read_reduced_state(const std::string & _file_base,
                                                    const bool _read_binary_data)
{
  // The reading mode: DECODE for binary, READ for ASCII
  XdrMODE _mode = _read_binary_data ? DECODE : READ;

  // The suffix to use for all the files that are written out
  const std::string _suffix = _read_binary_data ? ".xdr" : ".dat";

  // The string stream we'll use to make the file names
  std::ostringstream _file_name;

  // // First, find out how many basis functions we had when Greedy terminated
  // unsigned int n_bfs;

  _file_name << _file_base << _suffix;
  check_file_exists(_file_name.str());

  Xdr state_in(_file_name.str(), _mode);

  Number value;
  state_in >> value;
  _state_vector_reduced.push_back(value);
}

void
DwarfElephantReducedToFullState::check_file_exists(const std::string & _file_name)
{
  if (!std::ifstream(_file_name.c_str()))
    mooseError("The file: " + _file_name + " is not existing.");
}
