/**
 * This UserObject implements the Offline stage of the RB method.
 */

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantOfflineStage.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantOfflineStage>()
{
    InputParameters params = validParams<GeneralUserObject>();

    params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
    params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not.");
    params.addParam<bool>("compliant", true, "Determines whether F is equal to the output vector or not.");
    params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
    params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
    params.addParam<bool>("online_stage", false, "Determines whether the Online stage will be calculated or not.");
    params.addParam<std::string>("system","nl0","The name of the system that should be read in.");
    params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");
    params.addParam<Real>("mu_bar", 1., "Value for mu-bar");
//    params.addRequiredParam<unsigned int>("online_N","The number of basis functions that is used in the Reduced Basis solve during the Online Stage.");
    params.addRequiredParam<std::vector<Real>>("online_mu", "Current values of the different layers for which the RB Method is solved.");
    params.addRequiredParam<FunctionName>("cache_boundaries", "");

    return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
DwarfElephantOfflineStage::DwarfElephantOfflineStage(const InputParameters & params):
    GeneralUserObject(params),
    _use_displaced(getParam<bool>("use_displaced")),
    _store_basis_functions(getParam<bool>("store_basis_functions")),
    _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
    _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
    _compliant(getParam<bool>("compliant")),
    _online_stage(getParam<bool>("online_stage")),
    _system_name(getParam<std::string>("system")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
    _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject")),
    _function(&getFunction("cache_boundaries")),
    _mesh_ptr(&_fe_problem.mesh()),
    _subdomain_ids(_mesh_ptr->meshSubdomains()),
    _mu_bar(getParam<Real>("mu_bar")),
//    _online_N(getParam<unsigned int>("online_N")),
    _online_mu_parameters(getParam<std::vector<Real>>("online_mu"))
{
  _cache_boundaries = dynamic_cast<CacheBoundaries *>(_function);
}

void
DwarfElephantOfflineStage::setAffineMatrices()
{
//  if (_initialize_rb_system._qa > 1)
////  {
//   _initialize_rb_system._inner_product_matrix -> close();
    for(std::set<SubdomainID>::const_iterator it = _subdomain_ids.begin();
            it != _subdomain_ids.end(); it++)
    {
      _cache_boundaries->setCachedSubdomainStiffnessMatrixContributions(*_initialize_rb_system._jacobian_subdomain[*it], *it);
      _initialize_rb_system._jacobian_subdomain[*it] ->close();
//      _initialize_rb_system._inner_product_matrix->add(1., *_initialize_rb_system._jacobian_subdomain[*it]);
    }
//   }

//   else if (_initialize_rb_system._qa==1)
//   {
//     _initialize_rb_system._jacobian_subdomain[_initialize_rb_system._qa-1] -> close();
//     _initialize_rb_system._jacobian_subdomain[_initialize_rb_system._qa-1]->add(1, *_sys.matrix);
//   }
//
//    _cache_boundaries->setCachedStiffnessMatrixContributions(*_initialize_rb_system._inner_product_matrix);
    _initialize_rb_system._inner_product_matrix -> close();
    _initialize_rb_system._inner_product_matrix->add(-1., *_sys.matrix);
}

void
DwarfElephantOfflineStage::transferAffineVectors()
{
    // Transfer the vectors
    // Transfer the data for the F vectors.
    for(unsigned int _q=0; _q<_initialize_rb_system._qf; _q++)
    {
      _cache_boundaries->setCachedSubdomainResidual(*_initialize_rb_system._residuals[_q], _q);
      _initialize_rb_system._residuals[_q]->close();
    }

    // Transfer the data for the output vectors.
    if (_compliant)
    {
        for(unsigned int _q=0; _q<_initialize_rb_system._ql; _q++)
        {
          _cache_boundaries->setCachedSubdomainResidual(*_initialize_rb_system._outputs[_q], _q);
          _initialize_rb_system._outputs[_q]->close();
        }
    }
    else if (!_compliant)
        mooseError("Currently, the code handles the compliant case, only.");
}

void
DwarfElephantOfflineStage::offlineStage()
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

    _initialize_rb_system._rb_con_ptr->print_basis_function_orthogonality();
}

void
DwarfElephantOfflineStage::setOnlineParameters()
{
    for (unsigned int  _q = 0; _q != _online_mu_parameters.size(); _q++)
    {
        std::string  _mu_name = "mu_" + std::to_string(_q);
        _rb_online_mu.set_value(_mu_name, _online_mu_parameters[_q]);
    }
}

void
DwarfElephantOfflineStage::initialize()
{
}

void
DwarfElephantOfflineStage::execute()
{
    // Build the RBEvaluation object
    // Required for both the Offline and Online stage.
    DwarfElephantRBEvaluation  _rb_eval(_mesh_ptr->comm());

    // Pass a pointer of the RBEvaluation object to the
    // RBConstruction object
    _initialize_rb_system._rb_con_ptr->set_rb_evaluation(_rb_eval);

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

    if(_online_stage)
    {

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
      _rb_eval.rb_solve(_online_N);

//      for (unsigned int _q = 0; _q != _initialize_rb_system._ql; _q++)
//        _console << "Output " << std::to_string(_q) << ": value = " << _rb_eval.RB_outputs[_q]
//        << ", error bound = " << _rb_eval.RB_output_error_bounds[_q] << std::endl;

      _rb_eval.read_in_basis_functions(*_initialize_rb_system._rb_con_ptr);
      _initialize_rb_system._rb_con_ptr->load_rb_solution();
      ExodusII_IO(_mesh_ptr->getMesh()).write_equation_systems("RB_sol.e", _es);
      _initialize_rb_system._rb_con_ptr->load_basis_function(0);
      ExodusII_IO(_mesh_ptr->getMesh()).write_equation_systems("bf0.e", _es);
    }
}

void
DwarfElephantOfflineStage::finalize()
{
}
