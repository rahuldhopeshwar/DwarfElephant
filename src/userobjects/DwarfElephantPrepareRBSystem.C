 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantPrepareRBSystem.h"

template<>
InputParameters validParams<DwarfElephantPrepareRBSystem>()
{
  InputParameters params = validParams<NodalUserObject>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");

  return params;
}

DwarfElephantPrepareRBSystem::DwarfElephantPrepareRBSystem(const InputParameters & params):
  NodalUserObject(params),
//  _dof_map(_fe_problem.getNonlinearSystemBase().dofMap()),
  _use_displaced(getParam<bool>("use_displaced")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _Aq_qa(SparseMatrix<Number>::build(_es.comm()).release()),
  _sys(_es.get_system<TransientNonlinearImplicitSystem>("nl0"))
{
}

void
DwarfElephantPrepareRBSystem::initialize()
{
//  // Initialize the data structures for the stiffness matrix.
//  _dof_map.attach_matrix(*_Aq_qa);
//  _Aq_qa->init();
//  _Aq_qa->zero();
//
////  PetscMatrix<Number> * petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_Aq_qa);
////  MatSetOption(petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//
//  _Aq_qa->close();
//
////  _assembly.addJacobian(*_Aq_qa);
////  _fe_problem.addJacobian(*_Aq_qa,_tid);
//  _Aq_qa = _sys.matrix;
////  _Aq_qa->add(1,*_Aq_qa);
////  _Aq_qa->add(1,*_sys.matrix);
}

void
DwarfElephantPrepareRBSystem::execute()
{
//  _fe_problem.addJacobian(*_Aq_qa, _tid);
//  _console << _assembly.node()->id() << std::endl;
}

void
DwarfElephantPrepareRBSystem::threadJoin(const UserObject & y)
{
}

void
DwarfElephantPrepareRBSystem::finalize()
{
}
