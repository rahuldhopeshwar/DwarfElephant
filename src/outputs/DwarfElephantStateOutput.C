// MOOSE includes
#include "NonlinearSystemBase.h"
#include "FEProblem.h"

#include "DwarfElephantStateOutput.h"

template <>
InputParameters
validParams<DwarfElephantStateOutput>()
{
  InputParameters params = validParams<FileOutput>();
  params.addParam<std::string>("system","rb0","The name of the system that should be read in.");
  params.addParam<bool>("use_rb", false, "Defines whether the RB or FE method is used.");
  return params;
}

DwarfElephantStateOutput::DwarfElephantStateOutput(const InputParameters & parameters) :
    FileOutput(parameters),
    _use_rb(getParam<bool>("use_rb")),
    _system_name(getParam<std::string>("system"))
{
  std::ifstream input_file(filename());

  if (input_file)
  {
    remove(filename().c_str());
  }
}

void
DwarfElephantStateOutput::output(const ExecFlagType & /*type*/)
{
  Moose::perf_log.push("StateOutput()", "Output");

  if(processor_id() == 0)
  {
    std::filebuf _buffer;
    _buffer.open(filename(), std::ios::app);
    std::ostream state_file(&_buffer);

    if(!_use_rb)
    {
      _es_ptr->get_system(_system_name).solution->print(state_file);
    } else {
      mooseError("Currently under development.");
    }

    _buffer.close();
  }
  Moose::perf_log.pop("StateOutput()", "Output");
}



std::string
DwarfElephantStateOutput::filename()
{
  return _file_base;
}
