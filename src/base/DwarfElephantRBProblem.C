 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBProblem.h"

#include "Assembly.h"
#include "AuxiliarySystem.h"
#include "MooseEigenSystem.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<DwarfElephantRBProblem>()
{
  InputParameters params = validParams<FEProblemBase>();

  return params;
}

DwarfElephantRBProblem::DwarfElephantRBProblem(const InputParameters & params):
  FEProblemBase(params),
  _nl_sys((new DwarfElephantSystem(*this, "rb0")))

{
    _nl = _nl_sys;
    _aux = new AuxiliarySystem(*this, "aux0");
    
    //_assembly = _rb_assembly;
    
    newRBAssemblyArray(*_nl_sys);
 
    newAssemblyArray(*_nl_sys);
//    initNullSpaceVectors(parameters, *_nl_sys);

    _eq.parameters.set<DwarfElephantRBProblem *>("_fe_problem") = this;
    
    _console << "FE PID: " << processor_id() << std::endl;
}

DwarfElephantRBProblem::~DwarfElephantRBProblem()
{
  FEProblemBase::deleteAssemblyArray();

  delete _nl;

  delete _aux;
}

void
DwarfElephantRBProblem::setInputParametersFEProblem(InputParameters & parameters)
{
  // set _fe_problem
  FEProblemBase::setInputParametersFEProblem(parameters);
  // set _fe_problem
  parameters.set<DwarfElephantRBProblem *>("_fe_problem") = this;
}

void
DwarfElephantRBProblem::solve()
{
  Moose::perf_log.push("constructRB()", "Execution");

#ifdef LIBMESH_HAVE_PETSC
  Moose::PetscSupport::petscSetOptions(*this); // Make sure the PETSc options are setup for this app
#endif

  if (_solve)
    _nl_sys->solve();

  Moose::perf_log.pop("constructRB()", "Execution");
}

void
DwarfElephantRBProblem::newRBAssemblyArray(NonlinearSystemBase & nl)
{
  unsigned int n_threads = libMesh::n_threads();
  _rb_assembly.resize(n_threads);
  for (unsigned int i = 0; i < n_threads; i++)
    _rb_assembly[i] = new DwarfElephantRBAssembly(nl, i);
}

//void
//DwarfElephantRBProblem::newAssemblyArray(NonlinearSystemBase & nl)
//{
//  unsigned int n_threads = libMesh::n_threads();
//  _assembly.resize(n_threads);
// for (unsigned int i = 0; i < n_threads; i++)
//    _assembly[i] = new DwarfElephantRBAssembly(nl, i);
//  
//  _rb_assembly.resize(n_threads);
//  for (unsigned int i = 0; i < n_threads; i++)
//    _rb_assembly[i] = new DwarfElephantRBAssembly(nl, i);
//    
//    _console << "PID: " << processor_id() << std::endl; 
//}
