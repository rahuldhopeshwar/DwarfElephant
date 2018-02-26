#ifndef DWARFELEPHANTALLELEMENTALVARIABLEVALUES_H
#define DWARFELEPHANTALLELEMENTALVARIABLEVALUES_H

#include "GeneralVectorPostprocessor.h"

// Forward Declarations
class DwarfElephantAllElementalVariableValues;

template <>
InputParameters validParams<DwarfElephantAllElementalVariableValues>();

/**
 * Get all of the elements that are intersected by a line
 */
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
