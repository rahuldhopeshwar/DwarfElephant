// MOOSE includes
#include "NonlinearSystemBase.h"
#include "FEProblem.h"

#include "DwarfElephantDakotaOutput.h"

template <>
InputParameters
validParams<DwarfElephantDakotaOutput>()
{
  // InputParameters params = validParams<CSV>();
  InputParameters params = validParams<FileOutput>();
  params.addRequiredParam<std::vector<PostprocessorName>>("postprocessor", "Defines the name of the postprocessor(s) you want to use.");
  params.addParam<std::string>("delimiter", "   ", "Defines the delimiter.");
  // params.addParam<bool>("sort_columns", false, "Toggle the sorting of columns alphabetically.");

  // Options for aligning csv output with whitespace padding
  // params.addParam<bool>(
  //     "align",
  //     false,
  //     "Align the outputted csv data by padding the numbers with trailing whitespace");
  // params.addParam<unsigned int>("precision", 14, "Set the output precision");

  return params;
}

DwarfElephantDakotaOutput::DwarfElephantDakotaOutput(const InputParameters & parameters) :
    // CSV(parameters),
    FileOutput(parameters),
    // _align(getParam<bool>("align")),
    // _precision(getParam<unsigned int>("precision")),
    _delimiter(getParam<std::string>("delimiter")),
    // _sort_columns(getParam<bool>("sort_columns")),
    _postprocessor_name(getParam<std::vector<PostprocessorName>>("postprocessor"))
{
}

// void
// DwarfElephantDakotaOutput::initialSetup(){
//   CSV::initialSetup();
// }
//
// void
// DwarfElephantDakotaOutput::outputPostprocessors()
// {
//   CSV::outputPostprocessors();
//   _write_all_table = true;
// }
//
// void
// DwarfElephantDakotaOutput::outputVectorPostprocessors()
// {
//   CSV::outputVectorPostprocessors();
//   _write_vector_table = true;
// }
//
// void
// DwarfElephantDakotaOutput::output(const ExecFlagType & type)
// {
//   CSV::output(type);
// }

void
DwarfElephantDakotaOutput::output(const ExecFlagType & /*type*/)
{
  Moose::perf_log.push("DakotaOutput()", "Output");
  // This result file enables the use of MOOSE as a forward simulator within Dakota.
  // Which output parameters are printed to the result file can be controlled over the MOOSE input file.

//  if (type == EXEC_TIMESTEP_END)
//  {
  if(processor_id() == 0){
    std::ofstream dakota_file;
    dakota_file.open(filename(), std::ios::app);

    for(unsigned int i = 0; i < _postprocessor_name.size(); i++)
      if(i < _postprocessor_name.size()-1)
        dakota_file << _problem_ptr->getPostprocessorValue(_postprocessor_name[i]) << " f"<< _delimiter;
      else
        dakota_file << _problem_ptr->getPostprocessorValue(_postprocessor_name[i]) << " f"<< std::endl;
 // }
 // std::string deleteline = "0 f";
 // std::string line;
 //
 // std::ifstream dakota_file_rb;
 // dakota_file_rb.open(_file_path + _result_file_name + ".out", std::ios::app);
 //
 // while (std::getline(dakota_file_rb,line))
 // {
 //   line.replace(line.find(deleteline),deleteline.length(),"");
 //   dakota_file << line << std::endl;
 // }
    dakota_file.close();
  }
  Moose::perf_log.pop("DakotaOutput()", "Output");
}



std::string
DwarfElephantDakotaOutput::filename()
{
  return _file_base;
}
