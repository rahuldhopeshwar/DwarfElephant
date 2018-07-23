 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBExecutioner.h"
#include "DwarfElephantAppTypes.h"

template<>
InputParameters validParams<DwarfElephantRBExecutioner>()
{
  InputParameters params = validParams<Steady>();
    params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addParam<bool>("offline_stage", true, "Determines whether the Offline stage will be calculated or not.");
  return params;
}

DwarfElephantRBExecutioner::DwarfElephantRBExecutioner(const InputParameters & params):
  Steady(params),
  _simulation_type(getParam<std::string>("simulation_type")),
  _offline_stage(getParam<bool>("offline_stage"))
{
}

void
DwarfElephantRBExecutioner::init()
{
  std::cout << "Starting DwarfElephantRBExecutioner::init()" << std::endl;

  std::cout << "Starting DwarfElephantRBExecutioner::recoverApp" << std::endl;
  if (_app.isRecovering())
  {
    _console << "\nCannot recover steady solves!\nExiting...\n" << std::endl;
    return;
  }
  std::cout << "Done DwarfElephantRBExecutioner::recoverApp" << std::endl;

  std::cout << "Starting DwarfElephantRBExecutioner::checkIntegrity()" << std::endl;
  checkIntegrity();
  std::cout << "Done DwarfElephantRBExecutioner::checkIntegrity()" << std::endl;

  std::cout << "Starting DwarfElephantRBExecutioner::_problem.initialSetup()" << std::endl;
  _problem.initialSetup();
  std::cout << "Done DwarfElephantRBExecutioner::_problem.initialSetup()" << std::endl;

  _problem.execute(EXEC_EIM);

  std::cout << "Starting DwarfElephantRBExecutioner::_problem.outputStep(EXEC_INITIAL)" << std::endl;
  _problem.outputStep(EXEC_INITIAL);
  std::cout << "Done DwarfElephantRBExecutioner::_problem.outputStep(EXEC_INITIAL)" << std::endl;

  std::cout << "Done DwarfElephantRBExecutioner::init()" << std::endl;
}

void
DwarfElephantRBExecutioner::execute()
{
  std::cout << "Started DwarfElephantRBExecutioner execution" << std::endl;
  if (_app.isRecovering())
    return;

  preExecute();

  _problem.advanceState();

  // first step in any steady state solve is always 1 (preserving backwards compatibility)
  _time_step = 1;
  _time = _time_step;                 // need to keep _time in sync with _time_step to get correct output

// #ifdef LIBMESH_ENABLE_AMR
//
//   // Define the refinement loop
//   unsigned int steps = _problem.adaptivity().getSteps();
//   for (unsigned int r_step=0; r_step<=steps; r_step++)
//   {
// #endif //LIBMESH_ENABLE_AMR
    preSolve();
	std::cout << "RBExecutioner Pre-solve done" << std::endl;
    _problem.timestepSetup();
	std::cout << "RBExecutioner timeStepSetup done" << std::endl;
    _problem.execute(EXEC_TIMESTEP_BEGIN);
	std::cout << "RBExecutioner execute(EXEC_TIMESTEP_BEGIN) done" << std::endl;
    _problem.outputStep(EXEC_TIMESTEP_BEGIN);
	std::cout << "RBExecutioner outputStep(EXEC_TIMESTEP_BEGIN) done" << std::endl;

    // Update warehouse active objects
    _problem.updateActiveObjects();
	std::cout << "RBExecutioner updateActiveObjects() done" << std::endl;
    if (_offline_stage)
      _problem.solve();
  	std::cout << "RBExecutioner solve() done" << std::endl;
    postSolve();
		std::cout << "RBExecutioner post-solve done" << std::endl;

    _problem.onTimestepEnd();
	std::cout << "RBExecutioner onTimestepEnd() done" << std::endl;
    _problem.execute(EXEC_TIMESTEP_END);
	std::cout << "RBExecutioner execute(EXEC_TIMESTEP_END) done" << std::endl;

    if(_simulation_type == "steady")
    {
      _problem.computeIndicators();
	  	std::cout << "RBExecutioner computeIndicators() done" << std::endl;
      _problem.computeMarkers();
	std::cout << "RBExecutioner computeMarkers() done" << std::endl;
	
      _problem.execute(EXEC_CUSTOM);
	  	std::cout << "RBExecutioner execute(EXEC_CUSTOM) done" << std::endl;
      _problem.outputStep(EXEC_TIMESTEP_END);
	  	std::cout << "RBExecutioner outputStep(EXEC_TIMESTEP_END) done" << std::endl;
      _problem.outputStep(EXEC_CUSTOM);
	  	std::cout << "RBExecutioner outputStep(EXEC_CUSTOM) done" << std::endl;
    }

// #ifdef LIBMESH_ENABLE_AMR
//     if (r_step != steps)
//     {
//       _problem.adaptMesh();
//     }
//   }
// #endif

  postExecute();
}
