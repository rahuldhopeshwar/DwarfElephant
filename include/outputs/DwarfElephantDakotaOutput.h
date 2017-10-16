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
  std::string _system_name;
  std::string _variable_of_interest;
};

#endif /* DWARFELEPHANTDAKOTAOUTPUT_H */
