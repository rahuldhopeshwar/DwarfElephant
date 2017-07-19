/**
 * This UserObject implements the Offline and Online stage of the RB method.
 */

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantOfflineOnlineStageSteadyState.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantOfflineOnlineStageSteadyState>()
{
    InputParameters params = validParams<GeneralUserObject>();

    params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
    params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not.");
    params.addParam<bool>("compliant", true, "Determines whether F is equal to the output vector or not.");
    params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
    params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
    params.addParam<bool>("offline_stage", false, "Determines whether the Offline stage will be calculated or not.");
    params.addParam<bool>("online_stage", false, "Determines whether the Online stage will be calculated or not.");
    params.addParam<bool>("offline_error_bound", false, "Determines which error bound is used.");
    params.addParam<std::string>("system","rb0","The name of the system that should be read in.");
    params.addRequiredParam<std::string>("exodus_file_name","The name of the Exodus output file.");
    params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");
    params.addParam<Real>("mu_bar", 1., "Value for mu-bar");
    params.addRequiredParam<std::vector<Real>>("online_mu", "Current values of the different layers for which the RB Method is solved.");

    return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantOfflineOnlineStageSteadyState::DwarfElephantOfflineOnlineStageSteadyState(const InputParameters & params):
    GeneralUserObject(params),
    _use_displaced(getParam<bool>("use_displaced")),
    _store_basis_functions(getParam<bool>("store_basis_functions")),
    _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
    _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
    _compliant(getParam<bool>("compliant")),
    _offline_stage(getParam<bool>("offline_stage")),
    _online_stage(getParam<bool>("online_stage")),
    _offline_error_bound(getParam<bool>("offline_error_bound")),
    _system_name(getParam<std::string>("system")),
    _exodus_file_name(getParam<std::string>("exodus_file_name")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
    _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject")),
    _mesh_ptr(&_fe_problem.mesh()),
    _subdomain_ids(_mesh_ptr->meshSubdomains()),
    _mu_bar(getParam<Real>("mu_bar")),
    _online_mu_parameters(getParam<std::vector<Real>>("online_mu")),
    _rb_problem(cast_ptr<DwarfElephantRBProblem *>(&_fe_problem))
{
}

void
DwarfElephantOfflineOnlineStageSteadyState::setAffineMatrices()
{
  PARALLEL_TRY
  {
   _initialize_rb_system._inner_product_matrix -> close();
    for(unsigned int _q=0; _q<_initialize_rb_system._qa; _q++)
    {
      _rb_problem->rbAssembly(_q).setCachedStiffnessMatrixContributions(*_initialize_rb_system._jacobian_subdomain[_q]);
      _initialize_rb_system._jacobian_subdomain[_q] ->close();
      _initialize_rb_system._inner_product_matrix->add(_mu_bar, *_initialize_rb_system._jacobian_subdomain[_q]);
    }
  }
  PARALLEL_CATCH;
}

void
DwarfElephantOfflineOnlineStageSteadyState::transferAffineVectors()
{
  PARALLEL_TRY
  {
    // Transfer the vectors
    // Transfer the data for the F vectors.
    for(unsigned int _q=0; _q<_initialize_rb_system._qf; _q++)
    {
      _rb_problem->rbAssembly(_q).setCachedResidual(*_initialize_rb_system._residuals[_q]);
      _initialize_rb_system._residuals[_q]->close();
    }

    // Transfer the data for the output vectors.
    for(unsigned int i=0; i < _initialize_rb_system._n_outputs; i++)
    {
      for(unsigned int _q=0; _q < _initialize_rb_system._ql[i]; _q++)
      {
        _rb_problem->rbAssembly(_q).setCachedResidual(*_initialize_rb_system._outputs[i][_q]);
        _initialize_rb_system._outputs[i][_q]->close();
//      *_initialize_rb_system._outputs[i][_q] /= _mesh_ptr->nNodes();
//          _initialize_rb_system._outputs[_q]->set(100, 17.5);
      }
    }
  }
  PARALLEL_CATCH;
}

void
DwarfElephantOfflineOnlineStageSteadyState::offlineStage()
{
    _initialize_rb_system._rb_con_ptr->train_reduced_basis();
    #if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataSerialization::RBEvaluationSerialization _rb_eval_writer(_initialize_rb_system._rb_con_ptr->get_rb_evaluation());
      _rb_eval_writer.write_to_file("rb_eval.bin");
    #else
      // Write the offline data to file (xdr format).
      _initialize_rb_system._rb_con_ptr->get_rb_evaluation().legacy_write_offline_data_to_files();
    #endif

    // If desired, store the basis functions (xdr format).
    if (_store_basis_functions)
    {
        _initialize_rb_system._rb_con_ptr->get_rb_evaluation().write_out_basis_functions(*_initialize_rb_system._rb_con_ptr);
    }

//    _initialize_rb_system._rb_con_ptr->print_basis_function_orthogonality();
}

void
DwarfElephantOfflineOnlineStageSteadyState::setOnlineParameters()
{
    for (unsigned int  _q = 0; _q != _online_mu_parameters.size(); _q++)
    {
        std::string  _mu_name = "mu_" + std::to_string(_q);
        _rb_online_mu.set_value(_mu_name, _online_mu_parameters[_q]);
    }
}

void
DwarfElephantOfflineOnlineStageSteadyState::initialize()
{
}

void
DwarfElephantOfflineOnlineStageSteadyState::execute()
{

    // Build the RBEvaluation object
    // Required for both the Offline and Online stage.
    DwarfElephantRBEvaluationSteadyState _rb_eval(_mesh_ptr->comm() , _fe_problem);
    // Pass a pointer of the RBEvaluation object to the
    // RBConstruction object
    _initialize_rb_system._rb_con_ptr->set_rb_evaluation(_rb_eval);

    if (_offline_stage)
    {
       // Transfer the affine vectors to the RB system.
       if(_skip_vector_assembly_in_rb_system)
        transferAffineVectors();


      // Transfer the affine matrices to the RB system.
      if(_skip_matrix_assembly_in_rb_system)
        setAffineMatrices();

      // Perform the offline stage.
      _console << std::endl;
      offlineStage();
      _console << std::endl;
    }

    if(_online_stage)
    {
      Moose::perf_log.push("onlineStage()", "Execution");
      #if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataDeserialization::RBEvaluationDeserialization _rb_eval_reader(_rb_eval);
      #else
      _rb_eval.legacy_read_offline_data_from_files();
      #endif

      setOnlineParameters();
      _rb_eval.set_parameters(_rb_online_mu);

      _console << "---- Online Stage ----" << std::endl;
      _rb_eval.print_parameters();

      _online_N = _initialize_rb_system._rb_con_ptr->get_rb_evaluation().get_n_basis_functions();

      if(_offline_error_bound)
       _initialize_rb_system._rb_con_ptr->get_rb_evaluation().evaluate_RB_error_bound = false;

      _rb_eval.rb_solve(_online_N);

//      for (unsigned int i = 0; i != _initialize_rb_system._n_outputs; i++)
//        for (unsigned int _q = 0; _q != _initialize_rb_system._ql[i]; _q++)
//          _console << "Output " << std::to_string(i) << ": value = " << _rb_eval.RB_outputs[i]
//          << ", error bound = " << _rb_eval.RB_output_error_bounds[i] << std::endl;

      Moose::perf_log.push("write_exodus()", "Execution");

      std::string _systems_for_print[] = {"RBSystem"};
      const std::set<std::string>  _system_names_for_print (_systems_for_print, _systems_for_print+sizeof(_systems_for_print)/sizeof(_systems_for_print[0]));
      
      _rb_eval.read_in_basis_functions(*_initialize_rb_system._rb_con_ptr);
      _initialize_rb_system._rb_con_ptr->load_rb_solution();     
      
     ExodusII_IO(_mesh_ptr->getMesh()).write_equation_systems(_exodus_file_name + ".e", _es, &_system_names_for_print);


//      _initialize_rb_system._rb_con_ptr->load_basis_function(0);
//      ExodusII_IO(_mesh_ptr->getMesh()).write_equation_systems("bf0.e", _es);
    }
}

void
DwarfElephantOfflineOnlineStageSteadyState::finalize()
{
}
