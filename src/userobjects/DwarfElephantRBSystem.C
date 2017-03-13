 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBSystem.h"

template<>
InputParameters validParams<DwarfElephantRBSystem>()
{
  InputParameters params = validParams<NodalUserObject>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("offline_stage", true, "Determines whether the Offline stage will be calculated or not.");
  params.addParam<bool>("online_stage", true, "Determines whether the Online stage will be calculated or not.");
  params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not.");
  params.addParam<bool>("F_equal_to_output", true, "Determines whether F is equal to the output vector or not.");
  params.addRequiredParam<unsigned int>("online_N","The number of basis functions that is used in the Reduced Basis solve during the Online Stage.");
  params.addRequiredParam<Real>("online_mu", "Current values of the differnt layer for which the RB Method is solved.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the input file. Required for the libMesh functions");
  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");

  return params;
}

DwarfElephantRBSystem::DwarfElephantRBSystem(const InputParameters & params):
  NodalUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _skip_matrix_assembly(true),
  _skip_vector_assembly(true),
  _offline_stage(getParam<bool>("offline_stage")),
  _online_stage(getParam<bool>("online_stage")),
  _F_equal_to_output(getParam<bool>("F_equal_to_output")),
  _store_basis_functions(getParam<bool>("store_basis_functions")),
  _online_N(getParam<unsigned int>("online_N")),
  _online_mu(getParam<Real>("online_mu")),
  _parameters_filename(getParam<std::string>("parameters_filename")),
  _system_name(getParam<std::string>("system")),
  _block_ids(this->blockIDs()),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _mesh_ptr(&_fe_problem.mesh())
{
}

void
DwarfElephantRBSystem::transferAffineOperators()
{
  // Transfer the data for the F vectors.
  for(unsigned int _q=0; _q<_qf; _q++)
    _rb_con_ptr->get_Fq(_q)->operator=(*_sys.rhs);

  // Transfer the data for the output vectors.
  if (_F_equal_to_output)
  {
    for(unsigned int _q=0; _q<_ql; _q++)
      _rb_con_ptr->get_output_vector(0,_q)->operator=(*_sys.rhs);
  }
  else if (!_F_equal_to_output)
    mooseError("Currently, the code handles the case F is equal to the output vector, only.");

  // Transfer the A matrices.
  for (unsigned int _q = 0; _q < _qa; _q++)
  {
    SparseMatrix<Number> * _Aq_qa =_rb_con_ptr->get_Aq(_q);

    PetscMatrix<Number> * petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_Aq_qa);
    MatSetOption(petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    for (unsigned int i=0; i != _sys.n_dofs(); i++)
    {
      const Node & _node_ref_i = _mesh_ptr->nodeRef(i);
      const std::set<SubdomainID> & _node_block_ids_i = _mesh_ptr->getNodeBlockIds(_node_ref_i);

      for (std::set<SubdomainID>::const_iterator _it_i = _node_block_ids_i.begin();
           _it_i != _node_block_ids_i.end(); ++_it_i)
      {
        if(*_node_block_ids_i.find(*_it_i) == _q)
        {
          for (unsigned int j=0; j != _sys.n_dofs(); j++)
          {
            const Node & _node_ref_j = _mesh_ptr->nodeRef(j);
            const std::set<SubdomainID> & _node_block_ids_j = _mesh_ptr->getNodeBlockIds(_node_ref_j);

            for (std::set<SubdomainID>::const_iterator _it_j = _node_block_ids_j.begin();
             _it_j != _node_block_ids_j.end(); ++_it_j)
            {
              if(*_node_block_ids_i.find(*_it_j) == _q)
              {
                _Aq_qa->set(i, j, _sys.matrix->operator()(i,j));
              }
            }
          }
        }
      }
    }
    _Aq_qa->close();
  }
}

void
DwarfElephantRBSystem::offlineStage()
{
  // This method performs the offline stage of the RB problem.

  // Computation of the reduced basis space.
//  _rb_con_ptr->train_reduced_basis();
//
//  // Wrtite the offline data to file (xdr format).
//  _rb_con_ptr->get_rb_evaluation().legacy_write_offline_data_to_files();
//
//  // If desired, store the basis functions (xdr format).
//  if (_store_basis_functions)
//  {
//    _rb_con_ptr->get_rb_evaluation().write_out_basis_functions(*_rb_con_ptr);
//  }
//
//  _rb_con_ptr->print_basis_function_orthogonality();
}

void
DwarfElephantRBSystem::onlineStage()
{
  ////    _rb_eval.legacy_read_offline_data_from_files();
////    RBParameters _online_mu;
////
////    _online_mu.set_value("mu0", _online_mu0_parameters);
////    _rb_eval.set_parameters(_online_mu);
////    _rb_eval.rb_solve(_online_N);
////
////    _rb_eval.print_parameters();
}

void
DwarfElephantRBSystem::performRBSystem()
{
  // Build the RBEvaluation object.
  // Required for both the Offline and Online stage.
  DwarfElephantRBEvaluation _rb_eval(_mesh_ptr->comm());

  // Pass a pointer of the RBEvaluation object to the
  // RBConstruction object
  _rb_con_ptr->set_rb_evaluation(_rb_eval);

  if(_offline_stage && _online_stage)
  {
    offlineStage();
    onlineStage();
  }
//  else if(!_offline_stage && _online_stage)
//  {
//    onlineStage();
//  }
//
//  else if(_offline_stage && !_online_stage)
//  {
//    offlineStage();
//  }
//  else if(!_offline_stage && !_online_stage)
//  {
//    _console << "Both the offline and online stage are set to false, \
//                 meaning the RBSystem does nothing."
//  }
}


void
DwarfElephantRBSystem::initialize()
{
  // Define the parameter file for the libMesh functions.
  GetPot infile (_parameters_filename);

  // Add a new equation system for the RB construction.
  _rb_con_ptr = &_es.add_system<DwarfElephantRBConstruction> ("RBSystem");

  // Intialization of the added equation system
  _rb_con_ptr->init();

  // Build the RBEvaluation object
  // Required for both the Offline and Online stage.
  DwarfElephantRBEvaluation _rb_eval(_mesh_ptr->comm());

  // Pass a pointer of the RBEvaluation object to the
  // RBConstruction object
  _rb_con_ptr->set_rb_evaluation(_rb_eval);

  _qa = _rb_con_ptr->get_rb_theta_expansion().get_n_A_terms();
  _qf = _rb_con_ptr->get_rb_theta_expansion().get_n_F_terms();
  _ql = _rb_con_ptr->get_rb_theta_expansion().get_n_output_terms(0);

  if (_offline_stage)
  {
    // Get and process the necessary input parameters for the
    // offline stage
    _rb_con_ptr->process_parameters_file(_parameters_filename);

    // Print the system informations for the RBConstruction system.
    _rb_con_ptr->print_info();

    // Initialize the RB construction. Note, we skip the matrix and vector
    // assembly, since this is already done by MOOSE.
    _rb_con_ptr->initialize_rb_construction(_skip_matrix_assembly, _skip_vector_assembly);

//    prepareRBTraining();
  }
}

void
DwarfElephantRBSystem::execute()
{
}

void
DwarfElephantRBSystem::threadJoin(const UserObject & y)
{
}

void
DwarfElephantRBSystem::finalize()
{
 transferAffineOperators();
 performRBSystem();
}
