#include "DwarfElephantComputeEIMInnerProductMatrixSteadyState.h"
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantAppTypes.h"

//#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<DwarfElephantComputeEIMInnerProductMatrixSteadyState>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredParam<NonlinearVariableName>("variable","The name of the variable that this UserObject operates on");

  // General Parameters
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("offline_stage", true, "Determines whether the Offline stage will be calculated or not.");
  params.addParam<bool>("compliant", true, "Determines whether F is equal to the output vector or not.");
  params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
  params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");

  // Reduced Basis Parameters
  params.addRequiredParam<unsigned int>("n_training_samples_RB", "Defines the number of training samples used in the RB Greedy.");
  params.addParam<bool>("deterministic_training_RB", false, "Determines whether the RB training set is generated deterministically or randomly.");
  params.addParam<unsigned int>("training_parameters_random_seed_RB", -1, "Defines the random seed for the generation of the RB traning set.");
  params.addParam<bool>("quiet_mode_RB", true, "Determines the what is printed to the console.");
  params.addRequiredParam<unsigned int>("N_max_RB", "Defines the maximum number of RB basis functions.");
  params.addParam<Real>("rel_training_tolerance_RB", 1.0e-4, "Defines the relative training tolerance for the RB Greedy.");
  params.addParam<Real>("abs_training_tolerance_RB", 1.0e-12, "Defines the relative training tolerance for the RB Greedy.");
  params.addParam<bool>("normalize_rb_bound_in_greedy", false, "Determines whether the normalized RB bound is used in the Greedy or not.");
  params.addRequiredParam<std::vector<std::string>>("parameter_names_RB", "Parameter names for the RB method.");
  params.addParam<std::vector<Real>>("parameter_min_values_RB", "Defines the lower bound of the RB parameter range.");
  params.addParam<std::vector<Real>>("parameter_max_values_RB", "Defines the upper bound of the RB parameter range.");
  params.addParam<std::vector<std::string>>("discrete_parameter_names_RB", "Discrete parameter names for the RB method.");
  params.addParam<std::vector<Real>>("discrete_parameter_values_RB", "Defines the list of discrete RB parameter values.");
  params.addParam<std::string>("system","rb0","The name of the system that should be read in.");

  params.addRequiredParam<UserObjectName>("initialize_rb_userobject","Name of the userobject used to initialize the RBEIM system");
  ExecFlagEnum & exec = params.set<ExecFlagEnum>("execute_on");
  exec.addAvailableFlags(EXEC_EIM);
  params.setDocString("execute_on", exec.getDocString());

  return params;
}

DwarfElephantComputeEIMInnerProductMatrixSteadyState::DwarfElephantComputeEIMInnerProductMatrixSteadyState(const InputParameters & parameters) : 
  ElementUserObject(parameters),
  MooseVariableInterface<Real>(this, false, "variable"),
  _use_displaced(getParam<bool>("use_displaced")),
  _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _skip_vector_assembly_in_rb_system(getParam<bool>("skip_vector_assembly_in_rb_system")),
  _offline_stage(getParam<bool>("offline_stage")),
  _compliant(getParam<bool>("compliant")),
  _deterministic_training_RB(getParam<bool>("deterministic_training_RB")),
  _quiet_mode_RB(getParam<bool>("quiet_mode_RB")),
  _normalize_RB_bound_in_greedy(getParam<bool>("normalize_rb_bound_in_greedy")),
  _n_training_samples_RB(getParam<unsigned int>("n_training_samples_RB")),
  _training_parameters_random_seed_RB(getParam<unsigned int>("training_parameters_random_seed_RB")),
  _N_max_RB(getParam<unsigned int>("N_max_RB")),
  _system_name(getParam<std::string>("system")),
  _continuous_parameters_RB(getParam<std::vector<std::string>>("parameter_names_RB")),
  _discrete_parameters_RB(getParam<std::vector<std::string>>("discrete_parameter_names_RB")),
  _rel_training_tolerance_RB(getParam<Real>("rel_training_tolerance_RB")),
  _abs_training_tolerance_RB(getParam<Real>("abs_training_tolerance_RB")),
  _continuous_parameter_min_values_RB(getParam<std::vector<Real>>("parameter_min_values_RB")),
  _continuous_parameter_max_values_RB(getParam<std::vector<Real>>("parameter_max_values_RB")),
  _var(*mooseVariable()),
  _test(_var.phi()),
  _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initialize_rb_userobject"))
{
  std::cout << "DwarfElephantComputeEIMInnerProductMatrixSteadyState object created" << std::endl;
}

void
DwarfElephantComputeEIMInnerProductMatrixSteadyState::initialize()
{
  std::cout << "Starting DwarfElephantComputeEIMInnerProductMatrixSteadyStateInitialize" << std::endl;
  //_initialize_rb_system._inner_product_matrix_eim -> init();
}

void
DwarfElephantComputeEIMInnerProductMatrixSteadyState::execute()
{
  std::cout << "Starting DwarfElephantComputeEIMInnerProductMatrixSteadyState::execute" << std::endl;
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _test.size(); _j++)
      _local_ke(_i, _j) += computeIntegral(_i, _j);
  //std::cout << "_local_ke matrix calculated using _test" << "_test[" << 0 << "][" << 0 << "] = " << _test[0][0] << std::endl;
  //_local_ke.print();
    // If simulation is steady
  //const DwarfElephantInitializeRBEIMSystemSteadyState & _initialize_rbeim_system = getUserObject<DwarfElephantInitializeRBEIMSystemSteadyState>("initial_rbeim_userobjet");

  if (_initialize_rb_system._offline_stage)
      _initialize_rb_system._inner_product_matrix_eim -> add_matrix(_local_ke, _var.dofIndices());
  // make provision for modifying diagonal values, if required
}

Real
DwarfElephantComputeEIMInnerProductMatrixSteadyState::getValue()
{
  //gatherSum(_integral_value);
  return 0;//_integral_value;
}

void
DwarfElephantComputeEIMInnerProductMatrixSteadyState::threadJoin(const UserObject & y)
{
  // Is not executed in parallel runs
  const DwarfElephantComputeEIMInnerProductMatrixSteadyState & pps = static_cast<const DwarfElephantComputeEIMInnerProductMatrixSteadyState &>(y);
  std::cout << "Executing threadJoin()" << std::endl;
  //_integral_value += pps._integral_value;
}

Real
DwarfElephantComputeEIMInnerProductMatrixSteadyState::computeIntegral(unsigned int _i, unsigned int _j)
{
  Real sum = 0;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    sum += _JxW[_qp] * _coord[_qp] * _test[_i][_qp] * _test[_j][_qp];
  return sum;
}

void
DwarfElephantComputeEIMInnerProductMatrixSteadyState::processRBParameters()
{

  // End setting parameter values for the RB construction object
  // Set the random seed for the RNG. By default -1 is set, meaning that std::time is used as a seed for the RNG.
  _initialize_rb_system._rb_con_ptr->set_training_random_seed(_training_parameters_random_seed_RB);

  // Set quiet mode.
  _initialize_rb_system._rb_con_ptr->set_quiet_mode(_quiet_mode_RB);

  // Initialization of the RB parameters.
  _initialize_rb_system._rb_con_ptr->set_Nmax(_N_max_RB);

  _initialize_rb_system._rb_con_ptr->set_rel_training_tolerance(_rel_training_tolerance_RB);
  _initialize_rb_system._rb_con_ptr->set_abs_training_tolerance(_abs_training_tolerance_RB);

  _initialize_rb_system._rb_con_ptr->set_normalize_rb_bound_in_greedy(_normalize_RB_bound_in_greedy);

  RBParameters _mu_min_RB;
  RBParameters _mu_max_RB;

  for (unsigned int i=0; i<_continuous_parameters_RB.size(); i++)
  {
    _mu_min_RB.set_value(_continuous_parameters_RB[i], _continuous_parameter_min_values_RB[i]);
    _mu_max_RB.set_value(_continuous_parameters_RB[i], _continuous_parameter_max_values_RB[i]);
  }

  for (unsigned int i=0; i<_discrete_parameters_RB.size(); i++)
  {
    _discrete_parameter_values_RB[_discrete_parameters_RB[i]] = _discrete_parameter_values_in_RB;
  }

  std::map<std::string,bool> _log_scaling;
  RBParameters::const_iterator it     = _mu_min_RB.begin();
  RBParameters::const_iterator it_end = _mu_min_RB.end();
  for ( ; it != it_end; ++it)
    {
      std::string _param_name = it->first;

      // For now, just set all entries to false.
      // TODO: Implement a decent way to specify log-scaling true/false
      // in the input text file
      _log_scaling[_param_name] = false;
    }

    _initialize_rb_system._rb_con_ptr->initialize_parameters(_mu_min_RB, _mu_max_RB, _discrete_parameter_values_RB);

   _initialize_rb_system._rb_con_ptr->initialize_training_parameters(_initialize_rb_system._rb_con_ptr->get_parameters_min(),
                                               _initialize_rb_system._rb_con_ptr->get_parameters_max(),
                                               _n_training_samples_RB,
                                               _log_scaling,
                                               _deterministic_training_RB);
 
}

void DwarfElephantComputeEIMInnerProductMatrixSteadyState::finalize()
{
  _initialize_rb_system._inner_product_matrix_eim -> close();
  _initialize_rb_system._eim_con_ptr->train_reduced_basis();
  #if defined(LIBMESH_HAVE_CAPNPROTO)
    RBDataSerialization::RBEvaluationSerialization rb_eim_eval_writer(*(_initialize_rb_system._eim_eval_ptr));
    rb_eim_eval_writer.write_to_file("rb_eim_eval.bin");
  #else
    _initialize_rb_system._eim_con_ptr -> get_rb_evaluation().legacy_write_offline_data_to_files("eim_data");
  #endif
  processRBParameters();
  _initialize_rb_system._eim_eval_ptr -> initialize_eim_theta_objects();
  _initialize_rb_system._rb_eval_ptr -> get_rb_theta_expansion().attach_multiple_F_theta(_initialize_rb_system._eim_eval_ptr -> get_eim_theta_objects());
  _initialize_rb_system._eim_con_ptr -> initialize_eim_assembly_objects();
  _initialize_rb_system._rb_con_ptr -> get_rb_assembly_expansion().attach_multiple_F_assembly(_initialize_rb_system._eim_con_ptr -> get_eim_assembly_objects());
  _initialize_rb_system._rb_con_ptr -> print_info();
  
  _initialize_rb_system._rb_con_ptr -> initialize_rb_construction(_skip_matrix_assembly_in_rb_system, _skip_vector_assembly_in_rb_system);
  // Train reduced basis will be called after the kernel assembles the RB affine matrices and vectors
  
  _initialize_rb_system.AssignAffineMatricesAndVectors();
}
