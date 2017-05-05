 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystem.h"

template<>
InputParameters validParams<DwarfElephantInitializeRBSystem>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("offline_stage", true, "Determines whether the Offline stage will be calculated or not.");
  params.addParam<bool>("compliant", true, "Determines whether F is equal to the output vector or not.");
  params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
  params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the input file. Required for the libMesh functions");
  params.addRequiredParam<std::string>("rb_variable","Name of the variable for the RB method.");

  return params;
}

DwarfElephantInitializeRBSystem::DwarfElephantInitializeRBSystem(const InputParameters & params):
  GeneralUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _offline_stage(getParam<bool>("offline_stage")),
  _compliant(getParam<bool>("compliant")),
  _parameters_filename(getParam<std::string>("parameters_filename")),
  _rb_variable_name(getParam<std::string>("rb_variable")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _mesh_ptr(&_fe_problem.mesh()),
  _exec_flags(this->execFlags())
{
}


void
DwarfElephantInitializeRBSystem::initVariable()
{
//  unsigned int var_num = _rb_con_ptr->add_variable(_rb_variable_name, libMesh::FIRST);
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

   _inner_product_matrix = _rb_con_ptr->get_inner_product_matrix();
   PetscMatrix<Number> * _petsc_inner_matrix = dynamic_cast<PetscMatrix<Number>* > (_inner_product_matrix);
   MatSetOption(_petsc_inner_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//   _inner_product_matrix->close();

   for (unsigned int _q=0; _q < _qa; _q++)
   {
     _jacobian_subdomain[_q] = _rb_con_ptr->get_Aq(_q);

     // Eliminates error message for the initialization of new non-zero entries
     // For the future: change SparseMatrix pattern (increases efficency)
     PetscMatrix<Number> * _petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_jacobian_subdomain[_q]);
     MatSetOption(_petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//     _jacobian_subdomain[_q]->close();
    }

    for (unsigned int _q=0; _q < _qf; _q++)
      _residuals[_q] = _rb_con_ptr->get_Fq(_q);

    for (unsigned int _q=0; _q < _ql; _q++)
      _outputs[_q] = _rb_con_ptr->get_output_vector(_q,0);
}

void
DwarfElephantInitializeRBSystem::initialize()
{
  if (_exec_flags[0]==EXEC_INITIAL)
  {
    // Define the parameter file for the libMesh functions.
    GetPot infile (_parameters_filename);

    // Add a new equation system for the RB construction.
    _rb_con_ptr = &_es.add_system<DwarfElephantRBConstructionSteadyState> ("RBSystem");

    initVariable();
    // Intialization of the added equation system
    _rb_con_ptr->init();

     DwarfElephantRBEvaluationSteadyState _rb_eval(_mesh_ptr->comm()); //, _cache_boundaries, _fe_problem);
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
