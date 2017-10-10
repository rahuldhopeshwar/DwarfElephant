#ifndef DWARFELEPHANTRBOUTPUT_H
#define DWARFELEPHANTRBOUTPUT_H

// MOOSE includes
#include "FileOutput.h"
#include "MooseMesh.h"

// Forward declerations
class DwarfElephantRBOutput;

template <>
InputParameters validParams<DwarfElephantRBOutput>();

class DwarfElephantRBOutput : public FileOutput
{
public:
  DwarfElephantRBOutput(const InputParameters & parameters);

  virtual void output(const ExecFlagType & type) override;

protected:
  MooseMesh * _mesh_ptr;
};

#endif /* DWARFELEPHANTRBOUTPUT_H */
