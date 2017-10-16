// MOOSE includes
#include "NonlinearSystemBase.h"
#include "FEProblem.h"

#include "DwarfElephantDakotaOutput.h"

template <>
InputParameters
validParams<DwarfElephantDakotaOutput>()
{
  InputParameters params = validParams<FileOutput>();
  params.addParam<std::string>("system", "RBSystem", "Name of the used system for the solve.");
  params.addParam<std::string>("variable_of_interest", "", "Variable that should be exported to the output file.");

  return params;
}

DwarfElephantDakotaOutput::DwarfElephantDakotaOutput(const InputParameters & parameters) :
    FileOutput(parameters),
    _system_name(getParam<std::string>("system")),
    _variable_of_interest(getParam<std::string>("variable_of_interest"))
{
}


void
DwarfElephantDakotaOutput::output(const ExecFlagType & /*type*/)
{
  // This result file enables the use of MOOSE as a forward simulator within Dakota.
  // Which output parameters are printed to the result file can be controlled over the MOOSE input file.

//  NonlinearSystemBase & _nl_sys = _problem_ptr->getNonlinearSystemBase();

  std::ofstream dakota_file;
  dakota_file.open(_file_base + ".txt", std::ios::app);

//  if(_system_name == "rb0")
//  dakota_file << *_es_ptr->get_system<DwarfElephantRBConstructionSteadyState>(_system_name).solution << std::endl;

/// Retrieving matrices
//  dakota_file << *_es_ptr->get_system<DwarfElephantRBConstructionSteadyState>(_system_name).matrix << std::endl;

//  SparseMatrix<Number> & _matrix = *_es_ptr->get_system<DwarfElephantRBConstructionSteadyState>(_system_name).matrix;

//  _matrix.print_matlab("Matrix_Matlab_format");
}
