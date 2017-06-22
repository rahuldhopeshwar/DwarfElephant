#include "DwarfElephantSystem.h"

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
    // calculate the load vectors
    _computing_initial_residual = true;
    _fe_problem.computeResidual(_transient_sys, *_current_solution, *_transient_sys.rhs);
    _computing_initial_residual = false;
    _transient_sys.rhs->close();
  }

  // calculate the stiffness matrices
  _fe_problem.computeJacobian(_transient_sys, *_current_solution, *_transient_sys.matrix);
}
