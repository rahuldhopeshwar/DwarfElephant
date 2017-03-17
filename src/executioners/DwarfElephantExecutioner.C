 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantExecutioner.h"

template<>
InputParameters validParams<DwarfElephantExecutioner>()
{
  InputParameters params = validParams<Steady>();

  params += validParams<BlockRestrictable>();

  return params;
}

DwarfElephantExecutioner::DwarfElephantExecutioner(const InputParameters & params):
  Steady(params),
  BlockRestrictable(params)
{
}
