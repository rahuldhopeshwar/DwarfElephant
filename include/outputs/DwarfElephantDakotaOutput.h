#ifndef DWARFELEPHANTDAKOTAOUTPUT_H
#define DWARFELEPHANTDAKOTAOUTPUT_H

// MOOSE includes
#include "FileOutput.h"

// Forward declerations
class DwarfElephantDakotaOutput;

template <>
InputParameters validParams<DwarfElephantDakotaOutput>();

class DwarfElephantDakotaOutput : public FileOutput
{
public:
  DwarfElephantDakotaOutput(const InputParameters & parameters);

  virtual void output(const ExecFlagType & type) override;
};

#endif /* DWARFELEPHANTDAKOTAOUTPUT_H */
