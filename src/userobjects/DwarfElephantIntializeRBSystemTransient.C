 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemTransient.h"

template<>
InputParameters validParams<DwarfElephantInitializeRBSystemTransient>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("offline_stage", true, "Determines whether the Offline stage will be calculated or not.");
  params.addParam<bool>("compliant", true, "Determines whether F is equal to the output vector or not.");
  params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
  params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
  params.addParam<std::string>("system","rb0","The name of the system that should be read in.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the input file. Required for the libMesh functions");
  params.addRequiredParam<FunctionName>("cache_boundaries", "");

  return params;
}

DwarfElephantInitializeRBSystemTransient::DwarfElephantInitializeRBSystemTransient(const InputParameters & params):
  GeneralUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _offline_stage(getParam<bool>("offline_stage")),
  _compliant(getParam<bool>("compliant")),
  _system_name(getParam<std::string>("system")),
  _parameters_filename(getParam<std::string>("parameters_filename")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _mesh_ptr(&_fe_problem.mesh()),
  _sys(&_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _exec_flags(this->execFlags()),
  _function(&getFunction("cache_boundaries"))
{
  _cache_boundaries = dynamic_cast<CacheBoundaries *>(_function);
}


void
DwarfElephantInitializeRBSystemTransient::initializeOfflineStage()
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

   // Define size of all new parameters.
   _jacobian_subdomain.resize(_qa);
   _residuals.resize(_qf);
   _outputs.resize(_n_outputs);

   for (unsigned int i=0; i < _n_outputs; i++)
     _outputs[i].resize(_ql[i]);

    // Get the correct matrices from the RB System.

   // Eliminates error message for the initialization of new non-zero entries
   // For the future: change SparseMatrix pattern (increases efficency)
   _inner_product_matrix = _rb_con_ptr->get_inner_product_matrix();
   PetscMatrix<Number> * _petsc_inner_matrix = dynamic_cast<PetscMatrix<Number>* > (_inner_product_matrix);
   MatSetOption(_petsc_inner_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

   for (unsigned int _q=0; _q < _qa; _q++)
   {
     _jacobian_subdomain[_q] = _rb_con_ptr->get_Aq(_q);

     // Eliminates error message for the initialization of new non-zero entries
     // For the future: change SparseMatrix pattern (increases efficency)
     PetscMatrix<Number> * _petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_jacobian_subdomain[_q]);
     MatSetOption(_petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }

    // Get the correct vectors from the RB System.
    for (unsigned int _q=0; _q < _qf; _q++)
      _residuals[_q] = _rb_con_ptr->get_Fq(_q);

    for (unsigned int i=0; i < _n_outputs; i++)
      for (unsigned int _q=0; _q < _ql[i]; _q++)
        _outputs[i][_q] = _rb_con_ptr->get_output_vector(i,_q);
}

void
DwarfElephantInitializeRBSystemTransient::initialize()
{
  if (_exec_flags[0]==EXEC_INITIAL)
  {
    // Define the parameter file for the libMesh functions.
    GetPot infile (_parameters_filename);

    // Add a new equation system for the RB construction.
    _rb_con_ptr = &_es.add_system<DwarfElephantRBConstructionTransient> ("RBSystem");

    // Intialization of the added equation system
    _rb_con_ptr->init();
    _es.update();

    DwarfElephantRBEvaluationTransient _rb_eval(_mesh_ptr->comm(), _fe_problem);
    // Pass a pointer of the RBEvaluation object to the
    // RBConstruction object
    _rb_con_ptr->set_rb_evaluation(_rb_eval);

    TransientRBThetaExpansion & _trans_theta_expansion = cast_ref<TransientRBThetaExpansion &>(_rb_con_ptr->get_rb_theta_expansion());

    // Get number of attached parameters.
    _n_outputs = _rb_con_ptr->get_rb_theta_expansion().get_n_outputs();
    _ql.resize(_n_outputs);
    _qa = _rb_con_ptr->get_rb_theta_expansion().get_n_A_terms();
    _qm = _trans_theta_expansion.get_n_M_terms();
    _qf = _rb_con_ptr->get_rb_theta_expansion().get_n_F_terms();

    for(unsigned int i=0; i < _n_outputs; i++)
     _ql[i] = _rb_con_ptr->get_rb_theta_expansion().get_n_output_terms(i);

    // Initialize required matrices and vectors.
    if (_offline_stage)
    {
      initializeOfflineStage();
    }
  }
}

void
DwarfElephantInitializeRBSystemTransient::execute()
{
}

void
DwarfElephantInitializeRBSystemTransient::finalize()
{
}
