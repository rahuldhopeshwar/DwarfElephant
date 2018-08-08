#ifndef DWARFELEPHANTREVERSELIFTINGFUNCTIONANDDIMENSIONALIZE_H
#define DWARFELEPHANTREVERSELIFTINGFUNCTIONANDDIMENSIONALIZE_H

#include "NodalVectorPostprocessor.h"

// Forward Declarations
class DwarfElephantReverseLiftingFunctionAndDimensionalize;

template <>
InputParameters validParams<DwarfElephantReverseLiftingFunctionAndDimensionalize>();

/**
 * Get all of the elements that are intersected by a line
 */
class DwarfElephantReverseLiftingFunctionAndDimensionalize : public NodalVectorPostprocessor
{
public:
  DwarfElephantReverseLiftingFunctionAndDimensionalize(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin (const UserObject &/*uo*/) override {}

protected:
  Function & _lifting_function;
  NumericVector<Number> * _nodal_solution;
  VectorPostprocessorValue & _nodal_solution_original;
  std::string _system;
  Real _reference_value_variable;
  bool _dimensionalize;
  bool _scale_and_add;
  bool _reverse_lifting_function;
};

#endif
