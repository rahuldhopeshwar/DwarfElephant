// MOOSE includes
#include "NonlinearSystemBase.h"
#include "FEProblem.h"

#include "DwarfElephantDakotaOutput.h"

template <>
InputParameters
validParams<DwarfElephantDakotaOutput>()
{
  InputParameters params = validParams<FileOutput>();

  return params;
}

DwarfElephantDakotaOutput::DwarfElephantDakotaOutput(const InputParameters & parameters) :
    FileOutput(parameters) {}


void
DwarfElephantDakotaOutput::output(const ExecFlagType & /*type*/)
{
  NonlinearSystemBase & _nl_sys = _problem_ptr->getNonlinearSystemBase();

  std::ofstream dakota_file;
  dakota_file.open(_file_base + ".txt", std::ios::app);

  // Header of the file
  dakota_file << "Input file for Dakota generated within the DwarfElephant Application" << std::endl;
  dakota_file << std::endl;
  dakota_file << "****************************************************************************" << std::endl;
  dakota_file << "* Description: This input file enables the use of MOOSE as a forward       *" << std::endl;
  dakota_file << "*              simulator within Dakota. Which input parameters are defined *" << std::endl;
  dakota_file << "*              can be controlled over the MOOSE input file.                *" << std::endl;
  dakota_file << "****************************************************************************" << std::endl;
  dakota_file << std::endl;
}
