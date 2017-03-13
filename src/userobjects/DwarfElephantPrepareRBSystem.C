 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantPrepareRBSystem.h"

template<>
InputParameters validParams<DwarfElephantPrepareRBSystem>()
{
  InputParameters params = validParams<NodalUserObject>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");

  return params;
}

DwarfElephantPrepareRBSystem::DwarfElephantPrepareRBSystem(const InputParameters & params):
  NodalUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _system_name(getParam<std::string>("system")),
  _block_ids(this->blockIDs()),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _dof_map(_fe_problem.getNonlinearSystemBase().dofMap()),
  _mesh_ptr(&_fe_problem.mesh())
{
}

void
DwarfElephantPrepareRBSystem::initialize()
{
  // Rise an error in case the UserObject is defined for more than one block.
  if (_block_ids.size()>1)
    mooseError("The Userobject \"DwarfElephantPrepareRBSystem\" needs to be defined for each block separately."
               " You defined it for more than one block. Please change your specifications in the Inputfile.");

  // Initialize the data structures for the load vector.
  _Fq_a[*_block_ids.begin()] = NumericVector<Number>::build(_es.comm());
  _Fq_a[*_block_ids.begin()]->init(_sys.n_dofs(), _sys.n_local_dofs(), false, PARALLEL);

  // Initialize the data structures for the stiffness matrix.
  _Aq_a[*_block_ids.begin()] = SparseMatrix<Number>::build(_es.comm()).release();
  _dof_map.attach_matrix(*_Aq_a[*_block_ids.begin()]);
  _Aq_a[*_block_ids.begin()]->init();
  _Aq_a[*_block_ids.begin()]->zero();

  PetscMatrix<Number> * petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_Aq_a[*_block_ids.begin()]);
  MatSetOption(petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//  petsc_matrix->close();

//  for (std::set<SubdomainID>::const_iterator _q_a = _blocks.begin();
//      _q_a != _blocks.end(); ++_q_a)
//  {
//    _Fq_a[*_q_a] = NumericVector<Number>::build(_es.comm());
//    _Fq_a[*_q_a]->init(_rb_con_ptr->n_dofs(), _rb_con_ptr->n_local_dofs(), false, PARALLEL);
//
//    _Aq_a[*_q_a] = SparseMatrix<Number>::build(_es.comm());
//    _dof_map.attach_matrix(*_Aq_a[*_q_a]);
//    _Aq_a[*_q_a]->init();
//    _Aq_a[*_q_a]->zero();
//  }
}

void
DwarfElephantPrepareRBSystem::execute()
{
  // Split the residual vector into their subdomain contributions
  Number _value_vector = _sys.rhs->el(_current_node->id());
  _Fq_a[*_block_ids.begin()]->set(_current_node->id(), _value_vector);

  // Split the system matrix into their subdomain contributions
  for (unsigned int i= 0; i != _sys.n_dofs(); ++i)
  {
    const Node & _node_ref = _mesh_ptr->nodeRef(i);
    const std::set<SubdomainID> & _node_block_ids = _mesh_ptr->getNodeBlockIds(_node_ref);

    for (std::set<SubdomainID>::const_iterator _it = _node_block_ids.begin();
         _it != _node_block_ids.end(); ++_it)
      if(*_node_block_ids.find(*_it) == *_block_ids.begin())
        _Aq_a[*_block_ids.begin()]->set(_current_node->id(), i, _sys.matrix->operator()(_current_node->id(),i));
  }
}

void
DwarfElephantPrepareRBSystem::threadJoin(const UserObject & y)
{
}

void
DwarfElephantPrepareRBSystem::finalize()
{
  _Aq_a[*_block_ids.begin()]->close();
  _Aq_a[*_block_ids.begin()]->print();
  _console << _current_node << std::endl;
}
