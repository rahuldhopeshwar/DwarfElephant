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
    
  if (_solve)
    _nl_sys->update();

  Moose::perf_log.pop("constructRB()", "Execution");
}

void
DwarfElephantRBProblem::newRBAssemblyArray(NonlinearSystemBase & nl)
{
  unsigned int subdomains = mesh().meshSubdomains().size();
  _rb_assembly.resize(subdomains);
  for (unsigned int i = 0; i < subdomains; i++)
    _rb_assembly[i] = new DwarfElephantRBAssembly(nl, i);
}
