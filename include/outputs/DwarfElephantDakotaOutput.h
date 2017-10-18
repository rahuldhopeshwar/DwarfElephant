#ifndef DWARFELEPHANTDAKOTAOUTPUT_H
#define DWARFELEPHANTDAKOTAOUTPUT_H

// MOOSE includes
#include "FileOutput.h"
#include "DwarfElephantInitializeRBSystemSteadyState.h"

// Forward declerations
class DwarfElephantDakotaOutput;

template <>
InputParameters validParams<DwarfElephantDakotaOutput>();

class DwarfElephantDakotaOutput : public FileOutput
{
public:
  DwarfElephantDakotaOutput(const InputParameters & parameters);

  virtual void output(const ExecFlagType & type) override;

protected:
  std::string _result_file_name;
  std::string _file_path;

  PostprocessorName _postprocessor_name;
};

#endif /* DWARFELEPHANTDAKOTAOUTPUT_H */
