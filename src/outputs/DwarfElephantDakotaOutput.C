// MOOSE includes
#include "NonlinearSystemBase.h"
#include "FEProblem.h"

#include "DwarfElephantDakotaOutput.h"

template <>
InputParameters
validParams<DwarfElephantDakotaOutput>()
{
  InputParameters params = validParams<FileOutput>();
  params.addParam<std::string>("result_file_name", "results", "Defines the name of the result file.");
  params.addParam<std::string>("file_path", "" , "Defines the path where the file should be saved.");
  params.addRequiredParam<PostprocessorName>("postprocessor", "Defines the name of the postprocessor you want to use.");

  return params;
}

DwarfElephantDakotaOutput::DwarfElephantDakotaOutput(const InputParameters & parameters) :
    FileOutput(parameters),
    _result_file_name(getParam<std::string>("result_file_name")),
    _file_path(getParam<std::string>("file_path")),
    _postprocessor_name(getParam<PostprocessorName>("postprocessor"))
{
}


void
DwarfElephantDakotaOutput::output(const ExecFlagType & /*type*/)
{
  // This result file enables the use of MOOSE as a forward simulator within Dakota.
  // Which output parameters are printed to the result file can be controlled over the MOOSE input file.

//  if (type == EXEC_TIMESTEP_END)
//  {
    std::ofstream dakota_file;
    dakota_file.open(_file_path + _result_file_name + ".out", std::ios::app);
    dakota_file << _problem_ptr->getPostprocessorValue(_postprocessor_name) << " f"<< std::endl;
//  }
//  std::string deleteline = "0 f";
//  std::string line;
//
//  std::ifstream dakota_file_rb;
//  dakota_file_rb.open(_file_path + _result_file_name + ".out", std::ios::app);
//
//  while (std::getline(dakota_file_rb,line))
//  {
//    line.replace(line.find(deleteline),deleteline.length(),"");
//    dakota_file << line << std::endl;
//  }
}
