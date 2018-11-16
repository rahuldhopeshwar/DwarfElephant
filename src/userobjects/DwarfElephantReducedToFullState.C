///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantReducedToFullState.h"

registerMooseObject("DwarfElephantApp", DwarfElephantReducedToFullState);

template<>
InputParameters validParams<DwarfElephantReducedToFullState>()
{
  InputParameters params = validParams<GeneralUserObject>();
  return params;
}

DwarfElephantReducedToFullState::DwarfElephantReducedToFullState(const InputParameters & params):
  GeneralUserObject(params)
{
}

void
DwarfElephantReducedToFullState::initialize()
{
}

void
DwarfElephantReducedToFullState::execute()
{
  DwarfElephantRBProblem * _rb_problem = cast_ptr<DwarfElephantRBProblem *>(&_fe_problem);
  _rb_problem->setReducedInitialCondition();
}

void
DwarfElephantReducedToFullState::finalize()
{
  _fe_problem.computeIndicators();
  _fe_problem.computeMarkers();

  _fe_problem.execute(EXEC_CUSTOM);
  _fe_problem.outputStep(EXEC_TIMESTEP_END);
  _fe_problem.outputStep(EXEC_CUSTOM);
}
