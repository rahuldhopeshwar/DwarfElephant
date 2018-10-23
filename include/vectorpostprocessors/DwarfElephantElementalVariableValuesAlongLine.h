#ifndef DWARFELEPHANTELEMENTALVARIABLEVALUESALONGLINE_H
#define DWARFELEPHANTELEMENTALVARIABLEVALUESALONGLINE_H

#include "GeneralVectorPostprocessor.h"

// Forward Declarations
class DwarfElephantElementalVariableValuesAlongLine;

template <>
InputParameters validParams<DwarfElephantElementalVariableValuesAlongLine>();

///Get all of the elements that are intersected by a line
class DwarfElephantElementalVariableValuesAlongLine : public GeneralVectorPostprocessor
{
public:
  DwarfElephantElementalVariableValuesAlongLine(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;

protected:
  VectorPostprocessorValue & _elem_ids;
  VectorPostprocessorValue & _elem_x;
  VectorPostprocessorValue & _elem_y;
  VectorPostprocessorValue & _elem_z;
  VectorPostprocessorValue & _elem_values;
  const VariableName & _var_name;
  Point _start;
  Point _end;
};

#endif
