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
  NumericVector<Number> * _current_rb_solution = _sys.current_local_solution.get();
  _current_rb_solution->zero();

  if (_fe_problem.solverParams()._type != Moose::ST_LINEAR)
  {
    // set the boundaries for the FEM solutions
    _computing_initial_residual = true;
    // In case your are using a MOOSE version older than April 19th 2018 uncomment the following line
    // _fe_problem.computeResidual(_transient_sys, *_current_solution, *_transient_sys.rhs);
    _fe_problem.computeResidualSys(_transient_sys, *_current_rb_solution, *_transient_sys.rhs);

    _computing_initial_residual = false;
    _transient_sys.rhs->close();
  }

  // set the counters for the iterations to zero
  _current_l_its.clear();
  _current_nl_its = 0;

  // Initialize the solution vector
  setInitialSolution();
  // setRBInitialSolution();

  if(_use_finite_differenced_preconditioner)
    setupFiniteDifferencedPreconditioner();

//  system().solve();

// calculate the stiffness matrices
// In case your are using a MOOSE version older than April 19th 2018 uncomment the following line
// _fe_problem.computeJacobian(_transient_sys, *_current_solution, *_transient_sys.matrix);
_fe_problem.computeJacobianSys(_transient_sys, *_current_rb_solution, *_transient_sys.matrix);

//  DwarfElephantRBProblem & _rb_problem = cast_ref<DwarfElephantRBProblem &>(_fe_problem);

//  std::vector<std::string> _kernel_names = _rb_problem.getKernelNames();

//  for(unsigned int i = 0; i < _kernel_names.size(); i++)
//  {
//    std::shared_ptr<DwarfElephantRBKernel> _rb_kernel_ptr = std::static_pointer_cast<DwarfElephantRBKernel> (_kernels.getActiveObject(_kernel_names[i]));
//    _console << _kernel_names[i] << std::endl;
//    _rb_kernel_ptr->computeOutput();
//  }
}

void
DwarfElephantSystem::setRBInitialSolution()
{
  // deactiveAllMatrixTags();
  //
  // NumericVector<Number> & initial_solution(solution());
  // if (_predictor.get() && _predictor->shouldApply())
  // {
  //   _predictor->apply(initial_solution);
  //   _fe_problem.predictorCleanup(initial_solution);
  // }
  //
  // // do nodal BC
  // ConstBndNodeRange & bnd_nodes = *_mesh.getBoundaryNodeRange();
  // for (const auto & bnode : bnd_nodes)
  // {
  //   BoundaryID boundary_id = bnode->_bnd_id;
  //   Node * node = bnode->_node;
  //
  //   if (node->processor_id() == processor_id())
  //   {
  //     // reinit variables in nodes
  //     _fe_problem.reinitNodeFace(node, boundary_id, 0);
  //
  //     if (_preset_nodal_bcs.hasActiveBoundaryObjects(boundary_id))
  //     {
  //         const auto & preset_bcs = _preset_nodal_bcs.getActiveBoundaryObjects(boundary_id);
  //          for (const auto & preset_bc : preset_bcs)
  //            preset_bc->computeValue(initial_solution);
  //      }
  //    }
  //  }
  //
  //  _sys.solution->close();
  //  update();
  //
  //  // Set constraint slave values
  //  setConstraintSlaveValues(initial_solution, false);
  //
  //   if (_fe_problem.getDisplacedProblem())
  //      setConstraintSlaveValues(initial_solution, true);
  //   }
}
