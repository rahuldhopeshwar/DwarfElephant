#ifndef DWARFELEPHANTREVERSELIFTINGFUNCTION_H
#define DWARFELEPHANTREVERSELIFTINGFUNCTION_H

#include "NodalVectorPostprocessor.h"

// Forward Declarations
class DwarfElephantReverseLiftingFunction;

template <>
InputParameters validParams<DwarfElephantReverseLiftingFunction>();

/**
 * Get all of the elements that are intersected by a line
 */
class DwarfElephantReverseLiftingFunction : public NodalVectorPostprocessor
{
public:
  DwarfElephantReverseLiftingFunction(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin (const UserObject &/*uo*/) override {}

protected:
  Function & _lifting_function;
  NumericVector<Number> * _nodal_solution;
  VectorPostprocessorValue & _nodal_solution_original;
  std::string _system;
  // Real _scale;
};

#endif
