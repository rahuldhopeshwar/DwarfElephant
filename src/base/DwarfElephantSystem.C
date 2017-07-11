#include "DwarfElephantSystem.h"
#include "DwarfElephantRBProblem.h"

DwarfElephantSystem::DwarfElephantSystem(FEProblemBase & _fe_problem, const std::string & name) :
    NonlinearSystem(_fe_problem, name)
{
}

DwarfElephantSystem::~DwarfElephantSystem()
{
}

void
DwarfElephantSystem::solve()
{
  if (_fe_problem.solverParams()._type != Moose::ST_LINEAR)
  {
    // set the boundaries for the FEM solutions
    _computing_initial_residual = true;
    _fe_problem.computeResidual(_transient_sys, *_current_solution, *_transient_sys.rhs);
    //_console << *_transient_sys.rhs << std::endl;
    _computing_initial_residual = false;
    _transient_sys.rhs->close();
  }
  
  // set the counters for the iterations to zero
  _current_l_its.clear();
  _current_nl_its = 0;
  
  // Initialize the solution vector
  setInitialSolution();
  
  if(_use_finite_differenced_preconditioner)
    setupFiniteDifferencedPreconditioner();
    
  if (_time_integrator)
  {
    _console << "Time Integrator" << std::endl;
  }

  // calculate the stiffness matrices
  _fe_problem.computeJacobian(_transient_sys, *_current_solution, *_transient_sys.matrix);
  
  //_fe_problem.computeResidual(_transient_sys, *_current_solution, *_transient_sys.rhs);
}
