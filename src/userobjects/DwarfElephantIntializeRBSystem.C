 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystem.h"

template<>
InputParameters validParams<DwarfElephantInitializeRBSystem>()
{
  InputParameters params = validParams<GeneralUserObject>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("offline_stage", true, "Determines whether the Offline stage will be calculated or not.");
  params.addParam<bool>("online_stage", true, "Determines whether the Online stage will be calculated or not.");
  params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not.");
  params.addParam<bool>("F_equal_to_output", true, "Determines whether F is equal to the output vector or not.");
  params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
  params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
  params.addRequiredParam<unsigned int>("online_N","The number of basis functions that is used in the Reduced Basis solve during the Online Stage.");
  params.addRequiredParam<Real>("online_mu", "Current values of the differnt layer for which the RB Method is solved.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the input file. Required for the libMesh functions");
  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");

  return params;
}

DwarfElephantInitializeRBSystem::DwarfElephantInitializeRBSystem(const InputParameters & params):
  GeneralUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _offline_stage(getParam<bool>("offline_stage")),
  _online_stage(getParam<bool>("online_stage")),
  _F_equal_to_output(getParam<bool>("F_equal_to_output")),
  _store_basis_functions(getParam<bool>("store_basis_functions")),
  _online_N(getParam<unsigned int>("online_N")),
  _online_mu(getParam<Real>("online_mu")),
  _parameters_filename(getParam<std::string>("parameters_filename")),
  _system_name(getParam<std::string>("system")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _mesh_ptr(&_fe_problem.mesh())
{
}

void
DwarfElephantInitializeRBSystem::transferAffineOperators(bool _skip_matrix_assembly_in_rb_system, bool _skip_vector_assembly_in_rb_system)
{
  // Transfer the vectors
  if (_skip_vector_assembly_in_rb_system)
  {
    // Transfer the data for the F vectors.
    for(unsigned int _q=0; _q<_qf; _q++)
      _rb_con_ptr->get_Fq(_q)->operator=(_sys.get_vector("Re_non_time"));

    // Transfer the data for the output vectors.
    if (_F_equal_to_output)
    {
      for(unsigned int _q=0; _q<_ql; _q++)
        _rb_con_ptr->get_output_vector(0,_q)->operator=(_sys.get_vector("Re_non_time"));
    }
    else if (!_F_equal_to_output)
      mooseError("Currently, the code handles the compliant case, only.");
  }

  if (_skip_matrix_assembly_in_rb_system)
  {
    // The stiffness matrices are transfered in the RBKernel class.

    // Transfer the inner product matrix
    _rb_con_ptr->get_inner_product_matrix()->close();
    _rb_con_ptr->get_inner_product_matrix()->add(1,*_sys.matrix);
  }
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
DwarfElephantInitializeRBSystem::performRBSystem()
{
  // Build the RBEvaluation object.
  // Required for both the Offline and Online stage.
  DwarfElephantRBEvaluation _rb_eval(_mesh_ptr->comm());

  // Pass a pointer of the RBEvaluation object to the
  // RBConstruction object
  _rb_con_ptr->set_rb_evaluation(_rb_eval);

  if(_offline_stage && _online_stage)
  {
    //offlineStage();
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
DwarfElephantInitializeRBSystem::initialize()
{
  // Define the parameter file for the libMesh functions.
  GetPot infile (_parameters_filename);

  // Add a new equation system for the RB construction.
  _rb_con_ptr = &_es.add_system<DwarfElephantRBConstruction> ("RBSystem");
//  _rb_con_ptr = &_es.add_system<DwarfElephantRBConstruction> ("RBSystem");

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

    _rb_con_ptr->initialize_rb_construction(_skip_matrix_assembly_in_rb_system, _skip_vector_assembly_in_rb_system);
  }
}

void
DwarfElephantInitializeRBSystem::execute()
{
}

//void
//DwarfElephantInitializeRBSystem::threadJoin(const UserObject & y)
//{
//}


void
DwarfElephantInitializeRBSystem::finalize()
{
}
