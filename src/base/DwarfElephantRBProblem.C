 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBProblem.h"

template<>
InputParameters validParams<DwarfElephantRBProblem>()
{
  InputParameters params = validParams<FEProblem>();

  return params;
}

DwarfElephantRBProblem::DwarfElephantRBProblem(const InputParameters & params):
  FEProblem(params)
{
}
