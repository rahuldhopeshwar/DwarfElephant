 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantOfflineStage.h"

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
  params.addRequiredParam<std::string>("residual_name","Name of the residual vector that is retrieved. The name is either '_Re_non_time' or '_Re_time'.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");
  params.addRequiredParam<unsigned int>("online_N","The number of basis functions that is used in the Reduced Basis solve during the Online Stage.");
  params.addRequiredParam<std::vector<Real>>("online_mu", "Current values of the different layers for which the RB Method is solved.");
  params.addRequiredParam<AuxVariableName>("variable","Name of the variable for which the RB method will be performed.");

  return params;
}

DwarfElephantOfflineStage::DwarfElephantOfflineStage(const InputParameters & params):
  GeneralUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _store_basis_functions(getParam<bool>("store_basis_functions")),
  _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _compliant(getParam<bool>("compliant")),
  _online_stage(getParam<bool>("online_stage")),
  _system_name(getParam<std::string>("system")),
  _residual_name(getParam<std::string>("residual_name")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject")),
  _mesh_ptr(&_fe_problem.mesh()),
  _subdomain_ids(_mesh_ptr->meshSubdomains()),
  _online_N(getParam<unsigned int>("online_N")),
  _online_mu_parameters(getParam<std::vector<Real>>("online_mu")),
  _variable_name(params.get<AuxVariableName>("variable")),
  _variable(&_fe_problem.getVariable(_tid, _variable_name))
{
//  _variable = &_fe_problem.getVariable(_tid, _variable_name);
}

void
DwarfElephantOfflineStage::subdomainSetup()
{
  _console << "subdomain" << std::endl;
}

void
DwarfElephantOfflineStage::setInnerProductMatrix()
{
//  SparseMatrix<Number> * _inner_product_matrix = _initialize_rb_system._rb_con_ptr->get_inner_product_matrix();
//  _inner_product_matrix->close();

//  for(std::set<SubdomainID>::const_iterator it = _subdomain_ids.begin();
//      it != _subdomain_ids.end(); it++)
//  {
//     _initialize_rb_system._rb_con_ptr->get_Aq(*it)->close();
//     _inner_product_matrix->add(1,*_initialize_rb_system._rb_con_ptr->get_Aq(*it));
//  }

for(std::set<SubdomainID>::const_iterator it = _subdomain_ids.begin();
      it != _subdomain_ids.end(); it++)
  {
     _initialize_rb_system._rb_con_ptr->get_Aq(*it)->close();
  }
  _initialize_rb_system._inner_product_matrix -> close();
}

void
DwarfElephantOfflineStage::transferAffineVectors()
{
  // Transfer the vectors
  // Transfer the data for the F vectors.
  for(unsigned int _q=0; _q<_initialize_rb_system._qf; _q++)
    _initialize_rb_system._residuals[_q]->operator=(_sys.get_vector("RHS Vector"));

  // Transfer the data for the output vectors.
  if (_compliant)
  {
    for(unsigned int _q=0; _q<_initialize_rb_system._ql; _q++)
      _initialize_rb_system._outputs[_q]->operator=(_sys.get_vector("RHS Vector"));
  }
  else if (!_compliant)
    mooseError("Currently, the code handles the compliant case, only.");

//
//  if (_skip_matrix_assembly_in_rb_system)
//  {
//    // The stiffness matrices are transfered in the RBKernel class.
//
//    // Transfer the inner product matrix
//    _initialize_rb_system._rb_con_ptr->get_inner_product_matrix()->close();
//    _initialize_rb_system._rb_con_ptr->get_inner_product_matrix()->add(1,*_sys.matrix);
//  }
}

void
DwarfElephantOfflineStage::offlineStage()
{
  // This method performs the offline stage of the RB problem.

  // Computation of the reduced basis space.
  _initialize_rb_system._rb_con_ptr->train_reduced_basis();

  // Write the offline data to file (xdr format).
  _initialize_rb_system._rb_con_ptr->get_rb_evaluation().legacy_write_offline_data_to_files();

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
  if (_residual_name != "Re_non_time" && _residual_name != "Re_time")
    mooseError ("You have choosen an invalid residual_name. Valid names are 'Re_non_time' and 'Re_time'.");
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

  if(_skip_vector_assembly_in_rb_system)
    transferAffineVectors();

  if(_skip_matrix_assembly_in_rb_system)
    setInnerProductMatrix();

  _console << std::endl;
  offlineStage();
  _console << std::endl;

  if(_online_stage)
  {
//    _rb_eval.legacy_read_offline_data_from_files();
    setOnlineParameters();
    _rb_eval.set_parameters(_rb_online_mu);

    _console << "---- Online Stage ----" << std::endl;
    _rb_eval.print_parameters();

    _rb_eval.rb_solve(_online_N);

    for (unsigned int _q = 0; _q != _initialize_rb_system._ql; _q++)
      _console << "Output " << std::to_string(_q) << ": value = " << _rb_eval.RB_outputs[_q]
      << ", error bound = " << _rb_eval.RB_output_error_bounds[_q] << std::endl;

    _rb_eval.read_in_basis_functions(*_initialize_rb_system._rb_con_ptr);
    _console << *_initialize_rb_system._rb_con_ptr->solution << std::endl;
//    _fe_problem.getNonlinearSystemBase().computeJacobian(*_initialize_rb_system._jacobian_subdomain[0]);
//    _initialize_rb_system._jacobian_subdomain[0]->close();
    //_initialize_rb_system._jacobian_subdomain[0]->get_diagonal(_sys.get_vector(7));
  //  _console << _sys.get_vector(7) << std::endl;
//      _console << *_sys.rhs << std::endl;
//      _console << *_sys.solution << std::endl;
  }
}

void
DwarfElephantOfflineStage::finalize()
{
}
