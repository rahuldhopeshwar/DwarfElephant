 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystem.h"

template<>
InputParameters validParams<DwarfElephantInitializeRBSystem>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("offline_stage", true, "Determines whether the Offline stage will be calculated or not.");
//  params.addParam<bool>("online_stage", true, "Determines whether the Online stage will be calculated or not.");
//  params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not.");
  params.addParam<bool>("compliant", true, "Determines whether F is equal to the output vector or not.");
  params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
  params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
//  params.addRequiredParam<unsigned int>("online_N","The number of basis functions that is used in the Reduced Basis solve during the Online Stage.");
//  params.addRequiredParam<Real>("online_mu", "Current values of the differnt layer for which the RB Method is solved.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the input file. Required for the libMesh functions");
//  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");
  params.addRequiredParam<AuxVariableName>("variable","Name of the variable for which the RB method will be performed.");

  return params;
}

DwarfElephantInitializeRBSystem::DwarfElephantInitializeRBSystem(const InputParameters & params):
  GeneralUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _offline_stage(getParam<bool>("offline_stage")),
//  _online_stage(getParam<bool>("online_stage")),
  _compliant(getParam<bool>("compliant")),
//  _store_basis_functions(getParam<bool>("store_basis_functions")),
//  _online_N(getParam<unsigned int>("online_N")),
//  _online_mu(getParam<Real>("online_mu")),
  _parameters_filename(getParam<std::string>("parameters_filename")),
//  _system_name(getParam<std::string>("system")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
//  _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _mesh_ptr(&_fe_problem.mesh()),
  _exec_flags(this->execFlags()),
  _variable_name(params.get<AuxVariableName>("variable"))
//  _variable_name_lib(getParam<std::string>("variable"))
{
  _variable = &_fe_problem.getVariable(_tid,_variable_name);
}

void
DwarfElephantInitializeRBSystem::initializeOfflineStage()
{
  // Get and process the necessary input parameters for the
  // offline stage
  _rb_con_ptr->process_parameters_file(_parameters_filename);

  // Print the system informations for the RBConstruction system.
  _rb_con_ptr->print_info();

  // Initialize the RB construction. Note, we skip the matrix and vector
  // assembly, since this is already done by MOOSE.
  _rb_con_ptr->initialize_rb_construction(_skip_matrix_assembly_in_rb_system, _skip_vector_assembly_in_rb_system);

   // Save the A's, F's and output vectors from the RBConstruction class in pointers.
   // This additional saving of the pointers is required because otherwise a the RBEvaluation object has
   // to be set again in the RBKernel.
   _jacobian_subdomain.resize(_qa);
   _residuals.resize(_qf);
   _outputs.resize(_ql);

  for (unsigned int _q=0; _q < _qa; _q++)
  {
     _jacobian_subdomain[_q] = _rb_con_ptr->get_Aq(_q);

     // Eliminates error message for the initialization of new non-zero entries
     // For the future: change SparseMatrix pattern (increases efficency)
     PetscMatrix<Number> * _petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_jacobian_subdomain[_q]);
      MatSetOption(_petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   }

   for (unsigned int _q=0; _q < _ql; _q++)
     _residuals[_q] = _rb_con_ptr->get_Fq(_q);


   for (unsigned int _q=0; _q < _ql; _q++)
     _outputs[_q] = _rb_con_ptr->get_output_vector(0,_q);
}

void
DwarfElephantInitializeRBSystem::onlineStage()
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
DwarfElephantInitializeRBSystem::initialize()
{
  if (_exec_flags[0]==EXEC_INITIAL)
  {
    // Define the parameter file for the libMesh functions.
    GetPot infile (_parameters_filename);

    // Add a new equation system for the RB construction.
    _rb_con_ptr = &_es.add_system<DwarfElephantRBConstruction> ("RBSystem");

//    _rb_con_ptr->add_variable ("u", libMesh::FIRST);

    // Intialization of the added equation system
    _rb_con_ptr->init();

    // Build the RBEvaluation object
    // Required for both the Offline and Online stage.
    DwarfElephantRBEvaluation  _rb_eval(_mesh_ptr->comm());

    // Pass a pointer of the RBEvaluation object to the
    // RBConstruction object
    _rb_con_ptr->set_rb_evaluation(_rb_eval);

    _qa = _rb_con_ptr->get_rb_theta_expansion().get_n_A_terms();
    _qf = _rb_con_ptr->get_rb_theta_expansion().get_n_F_terms();
    _ql = _rb_con_ptr->get_rb_theta_expansion().get_n_output_terms(0);


    if (_offline_stage)
    {
      initializeOfflineStage();
    }
  }
}

void
DwarfElephantInitializeRBSystem::execute()
{
}

void
DwarfElephantInitializeRBSystem::finalize()
{
}
