// MOOSE includes
#include "NonlinearSystemBase.h"
#include "FEProblem.h"

#include "DwarfElephantDakotaOutput.h"

template <>
InputParameters
validParams<DwarfElephantDakotaOutput>()
{
  InputParameters params = validParams<FileOutput>();
  // params.addRequiredParam<std::string>("postprocessor", "Defines the name of the postprocessor(s) you want to use.");
  params.addRequiredParam<std::vector<PostprocessorName>>("postprocessor", "Defines the name of the postprocessor(s) you want to use.");
  params.addParam<std::string>("delimiter", "   ", "Defines the delimiter.");
  return params;
}

DwarfElephantDakotaOutput::DwarfElephantDakotaOutput(const InputParameters & parameters) :
    FileOutput(parameters),
    _postprocessor_name(getParam<std::vector<PostprocessorName>>("postprocessor")),
    _delimiter(getParam<std::string>("delimiter"))
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
    dakota_file.open(filename() + ".out", std::ios::app);

    for(unsigned int i = 0; i < _postprocessor_name.size(); i++)
      if(i < _postprocessor_name.size()-1)
        dakota_file << _problem_ptr->getPostprocessorValue(_postprocessor_name[i]) << " f"<< _delimiter;
      else
        dakota_file << _problem_ptr->getPostprocessorValue(_postprocessor_name[i]) << " f"<< std::endl;
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
