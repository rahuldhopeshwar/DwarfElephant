#ifndef DWARFELEPHANTALLELEMENTALVARIABLEVALUES_H
#define DWARFELEPHANTALLELEMENTALVARIABLEVALUES_H

#include "GeneralVectorPostprocessor.h"

// Forward Declarations
class DwarfElephantAllElementalVariableValues;

template <>
InputParameters validParams<DwarfElephantAllElementalVariableValues>();

/// This VectorPostprocessor returns all of the elemental variable values of the mesh.
class DwarfElephantAllElementalVariableValues : public GeneralVectorPostprocessor
{
public:
  DwarfElephantAllElementalVariableValues(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  // virtual void threadJoin(const UserObject & s) override;

protected:
  VectorPostprocessorValue & _elem_values;
  const VariableName & _var_name;
};

#endif
