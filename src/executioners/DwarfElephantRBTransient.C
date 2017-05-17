 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBTransient.h"

template<>
InputParameters validParams<DwarfElephantRBTransient>()
{
  InputParameters params = validParams<Transient>();

  return params;
}

DwarfElephantRBTransient::DwarfElephantRBTransient(const InputParameters & params):
  Transient(params)
{
}

void
DwarfElephantRBTransient::execute()
{
//  preExecute();
//
//  // NOTE: if you remove this line, you will see a subset of tests failing. Those tests might have a wrong answer and might need to be regolded.
//  // The reason is that we actually move the solution back in time before we actually start solving (which I think is wrong).  So this call here
//  // is to maintain backward compatibility and so that MOOSE is giving the same answer.  However, we might remove this call and regold the test
//  // in the future eventually.
//  if (!_app.isRecovering())
//    _problem.advanceState();
//
//  // Start time loop...
//  while (true)
//  {
//    if (_first != true)
//      incrementStepOrReject();
//
//    _first = false;
//
//    if (!keepGoing())
//      break;
//
//    preStep();
//    computeDT();
//    takeStep();
//    endStep();
//    postStep();
//
//    _steps_taken++;
//  }
//
//  if (!_app.halfTransient())
//    _problem.outputStep(EXEC_FINAL);
//  postExecute();
}
