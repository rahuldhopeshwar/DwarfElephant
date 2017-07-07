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
    // calculate the load vectors
    _computing_initial_residual = true;
    _fe_problem.computeResidual(_transient_sys, *_current_solution, *_transient_sys.rhs);
    //_console << *_transient_sys.rhs << std::endl;
    _computing_initial_residual = false;
    _transient_sys.rhs->close();
  }

  // calculate the stiffness matrices
  _fe_problem.computeJacobian(_transient_sys, *_current_solution, *_transient_sys.matrix);
  //_console << *_transient_sys.matrix << std::endl;

//  DwarfElephantRBConstructionSteadyState * _rb_sys = &_fe_problem.es().get_system<DwarfElephantRBConstructionSteadyState>("RBSystem");
//  DwarfElephantRBEvaluationSteadyState _rb_eval(_fe_problem.mesh().comm(),_fe_problem);
//  _rb_sys->set_rb_evaluation(_rb_eval);
//  _rb_sys->print_info();
//
//  DwarfElephantRBProblem * _rb_problem = cast_ptr<DwarfElephantRBProblem *>(&_fe_problem);
//
//    PARALLEL_TRY
//  {
//      //_rb_problem->rbAssembly(0).setCachedSubdomainStiffnessMatrixContributions(*_rb_sys->get_Aq(0),0);
////      _rb_sys->get_Aq(0)->close();
//  }
//  PARALLEL_CATCH;
}
