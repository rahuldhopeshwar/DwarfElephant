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
  virtual void threadJoin (const UserObject &/*uo*/) override;

protected:
  NumericVector<Number> * _nodal_solution;
  // std::unique_ptr<NumericVector<Number>> _nodal_solution_add_on;
  VectorPostprocessorValue & _nodal_solution_original;
  std::string _system;
  Real _reference_value_variable;
  Real _scale_lifting_function;
  Real _range;
  bool _dimensionalize;
  bool _scale_and_add;
  bool _kriging;
  bool _reverse_lifting_function;
  Function * _lifting_function_1;
  Function * _lifting_function_2;
};

#endif
