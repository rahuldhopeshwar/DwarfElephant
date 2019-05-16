#include "DwarfElephantRBTransient.h"

// MOOSE includes
#include "Factory.h"
#include "SubProblem.h"
#include "TimeStepper.h"
#include "MooseApp.h"
#include "Conversion.h"
#include "FEProblem.h"
#include "NonlinearSystem.h"
#include "Control.h"
#include "TimePeriod.h"
#include "MooseMesh.h"
#include "AllLocalDofIndicesThread.h"
#include "TimeIntegrator.h"

#include "libmesh/implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/numeric_vector.h"

// C++ Includes
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

registerMooseObject("DwarfElephantApp", DwarfElephantRBTransient);

template <>
InputParameters
validParams<DwarfElephantRBTransient>()
{
  InputParameters params = validParams<Transient>();
  return params;
}

DwarfElephantRBTransient::DwarfElephantRBTransient(const InputParameters & parameters)
  : Transient(parameters)
{
}

void
DwarfElephantRBTransient::execute()
{
  if (_app.isRecovering())
      return;

    preExecute();

    _problem.advanceState();

    // first step in any steady state solve is always 1 (preserving backwards compatibility)
    _t_step = 1;
    _time = _t_step;                 // need to keep _time in sync with _time_step to get correct output

  // #ifdef LIBMESH_ENABLE_AMR
  //
  //   // Define the refinement loop
  //   unsigned int steps = _problem.adaptivity().getSteps();
  //   for (unsigned int r_step=0; r_step<=steps; r_step++)
  //   {
  // #endif //LIBMESH_ENABLE_AMR
      preSolve();
      _problem.timestepSetup();
      _problem.execute(EXEC_TIMESTEP_BEGIN);
      _problem.outputStep(EXEC_TIMESTEP_BEGIN);

      // Update warehouse active objects
      _problem.updateActiveObjects();

      // if (_offline_stage)
        _problem.solve();
      postSolve();

      _problem.onTimestepEnd();
      // _problem.execute(EXEC_CUSTOM);
      _problem.execute(EXEC_TIMESTEP_END);
      // _problem.outputStep(EXEC_CUSTOM);

      // if(_simulation_type == "steady")
      // {
      //   // _problem.computeIndicators();
      //   // _problem.computeMarkers();
      //   //
      //   // _problem.execute(EXEC_CUSTOM);
      //   // _problem.outputStep(EXEC_TIMESTEP_END);
      //   // _problem.outputStep(EXEC_CUSTOM);
      // }

  // #ifdef LIBMESH_ENABLE_AMR
  //     if (r_step != steps)
  //     {
  //       _problem.adaptMesh();
  //     }
  //   }
  // #endif

      {
        TIME_SECTION(_final_timer)
        _problem.execute(EXEC_FINAL);
        _problem.outputStep(EXEC_FINAL);
      }

    postExecute();
}
