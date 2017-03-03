///---------------------------------INCLUDE---------------------------------
// Moose includes
#include "FEProblem.h"
#include "FileMesh.h"
#include "MooseApp.h"

// Moose includes (DwarfElephant package)
#include "RBOutput.h"

//#include "libmesh/rb_data_serialization.h"
///-------------------------------------------------------------------------
template<>

///----------------------------INPUT PARAMETERS-----------------------------
InputParameters validParams<RBOutput>()
{

  // Get the parameters from the parent object
  InputParameters params = validParams<AdvancedOutput< FileOutput > >();
  params += AdvancedOutput< FileOutput >::enableOutputTypes("scalar postprocessor input system_information");

  // RB Parameters
  params.addRequiredParam<bool>("offline_stage", "Determines whether the Offline stage will be calculated or not.");
  params.addRequiredParam<bool>("online_stage", "Determines whether the Online stage will be calculated or not.");
  params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not");
  params.addRequiredParam<unsigned int>("online_N","The number of basis functions that is used in the Reduced Basis solve during the Online Stage.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the input file. Required for the libMesh functions");
  params.addRequiredParam<Real>("online_mu", "Current values of the differnt layer for which the RB Method is solved.");
  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");

  //params.set<MultiMooseEnum>("execute_on", /*quiet_mode=*/true).push_back("initial timestep_begin linear nonlinear failed");
  params.set<MultiMooseEnum>("execute_system_information_on") = "initial timestep_end";
  params.set<MultiMooseEnum>("execute_postprocessors_on") = "initial timestep_end";
  params.set<MultiMooseEnum>("execute_scalars_on") = "initial timestep_end";

  params.set<std::string>("built_by_action") = "add_user_object";


  return params;
}

///-------------------------------------------------------------------------
RBOutput::RBOutput(const InputParameters & parameters) :
    AdvancedOutput<FileOutput>(parameters),
    parameters_filename(getParam<std::string>("parameters_filename")),
    _offline_stage(getParam<bool>("offline_stage")),
    _online_stage(getParam<bool>("online_stage")),
    _store_basis_functions(getParam<bool>("store_basis_functions")),
    _online_N(getParam<unsigned int>("online_N")),
    _online_mu(getParam<Real>("online_mu")),
    _skip_matrix_assembly(false),
    _skip_vector_assembly(false),
    _mesh_ptr(&_problem_ptr->mesh()),
    _system_name(getParam<std::string>("system")),
    _system(&_es_ptr->get_system(_system_name))
{
}

///-------------------------------------------------------------------------
void
RBOutput::prepareDataStructuresRB()
{
//  _console << std::endl;
//  unsigned int _n_es = _es_ptr->n_systems();
//  _console << "Number of EquationsSystems: " << _n_es << std::endl;
//
//  unsigned int _n_vars_es = _es_ptr->n_vars();
//  _console << "The equations system has " << _n_vars_es << " variables." << std::endl;
//
//  bool _has_desired_system = _es_ptr->has_system(_system_name);
//  _console << "The " << _system_name << " system exists? " << _has_desired_system << std::endl;
//  _console << std::endl;
//
//  _system = &_es_ptr->get_system(_system_name);
//
//  const std::string _name_sys = _system->name();
//  std::string _type_sys = _system->system_type();
//  unsigned int _number_sys = _system->number();
//  _console << "The name of the obtained system is " << _name_sys << ", the type is "
//           << _type_sys << ", and it has the number " << _number_sys << "." << std::endl;
//
//  unsigned int _n_vectors_sys = _system->n_vectors();
//  unsigned int _n_matrices_sys = _system->n_matrices();
//  unsigned int _n_vars_sys = _system->n_vars();
//  _console << "The " << _name_sys << " system has " << _n_vars_sys << " variables, "
//           << _n_vectors_sys << " vectors, and " << _n_matrices_sys << " matrices." << std::endl;
//
  const std::string & _vec0_sys = _system->vector_name(0);
  const std::string & _vec1_sys = _system->vector_name(1);
  const std::string & _vec2_sys = _system->vector_name(2);
  const std::string & _vec3_sys = _system->vector_name(3);
  const std::string & _vec4_sys = _system->vector_name(4);
  const std::string & _vec5_sys = _system->vector_name(5);
  const std::string & _vec6_sys = _system->vector_name(6);
  const std::string & _vec7_sys = _system->vector_name(7);
  _console << "Name of the first vector: " << _vec0_sys << std::endl;
  _console << "Name of the second vector: " << _vec1_sys << std::endl;
  _console << "Name of the third vector: " << _vec2_sys << std::endl;
  _console << "Name of the fourth vector: " << _vec3_sys << std::endl;
  _console << "Name of the fifth vector: " << _vec4_sys << std::endl;
  _console << "Name of the sixth vector: " << _vec5_sys << std::endl;
  _console << "Name of the seventh vector: " << _vec6_sys << std::endl;
  _console << "Name of the eigth vector: " << _vec7_sys << std::endl;

//  NumericVector<Number> & _residual_non_time = _system->get_vector(1);
//  std::string _residual_info = _residual_non_time.get_info();
//  _console << _residual_info << std::endl;

//  NonlinearImplicitSystem * _imp_sys = &_es_ptr->get_system<NonlinearImplicitSystem>(_system_name);
//  FEProblemBase * _p_ptr = _imp_sys->get_equation_systems().parameters.get<FEProblemBase *>("_fe_problem_base");
//  NonlinearSystem * _non_sys_ptr = &_p_ptr->getNonlinearSystem();
//  std::string _name_sys = _non_sys_ptr->name();
//  _console << _name_sys << std::endl;
//  SparseMatrix<Number> & _jacobian_non_time = _imp_sys->get_matrix("Preconditioner");
  ///-------------------------------------------------------------------------------------------------------

}

void
RBOutput::performRBSystemOld()
{

//  GetPot infile (parameters_filename);

  DwarfElephantRBConstruction & _rb_con =
    _es_ptr->add_system<DwarfElephantRBConstruction> ("RBSystem");

  // Intialization of the added equation system
  _es_ptr->init();
//
  // Build the RBEvaluation object
  // Required for both the Offline and Online stage.
  DwarfElephantRBEvaluation _rb_eval(_mesh_ptr->comm());

  // Pass a pointer of the RBEvaluation object to the
  // RBConstruction object
  _rb_con.set_rb_evaluation(_rb_eval);


  if(_offline_stage && _online_stage)
  {
    /// Offline stage
    // Initialize the RB Construction process.
    // In this method the required data structures are set up and performed.
    // Also the intial assembly of the "truth" affine expansion of the PDE
    // is implemented here.

    _rb_con.process_parameters_file(parameters_filename);

    _rb_con.print_info();

    _rb_con.initialize_rb_construction();

//    _rb_con.solve_for_matrix_and_rhs(*_rb_con.get_linear_solver(), _rb_con.get_matrix_for_output_dual_solves(), *_rb_con.rhs);

    // Computation of the reduced basis space by taking snapshots.
    _rb_con.train_reduced_basis();
////
//
//    NumericVector<Number> & _Fq = *_rb_con.get_Fq(0);
//    Real _min_Fq = _Fq.min();
//    Real _max_Fq = _Fq.max();
//
//    _console << "libMesh min: " << _min_Fq << std::endl;
//    _console << "libMesh max: " << _max_Fq << std::endl;
//
//
    _rb_con.get_rb_evaluation().legacy_write_offline_data_to_files();
//
    if (_store_basis_functions)
    {
      _rb_con.get_rb_evaluation().write_out_basis_functions(_rb_con);
    }

    _rb_con.print_basis_function_orthogonality();
//
//
////    /// Online stage
////    RBParameters _online_mu;
////
////    _online_mu.set_value("mu0", _online_mu0_parameters);
////    _rb_eval.set_parameters(_online_mu);
////    _rb_eval.rb_solve(_online_N);
////
////    _rb_eval.print_parameters();
  }
//////    elseif(_offline_stage==false && _online_stage==true)
//////     _rb_eval.legacy_read_offline_data_from_files();
//////    elseif(_offline_stage==true && _online_stage==false)
//////    elseif(_offline_stage==false && _online_stage==false
////    //return _rb_con;

  }

void
RBOutput::output(const ExecFlagType & type)
{
  if (type==EXEC_TIMESTEP_END)
  {
    performRBSystemOld();
  }
}
