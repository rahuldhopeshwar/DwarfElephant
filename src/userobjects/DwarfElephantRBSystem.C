 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBSystem.h"

template<>
InputParameters validParams<DwarfElephantRBSystem>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params += validParams<BlockRestrictable>();
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("offline_stage", true, "Determines whether the Offline stage will be calculated or not.");
  params.addParam<bool>("online_stage", true, "Determines whether the Online stage will be calculated or not.");
  params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not");
  params.addRequiredParam<unsigned int>("online_N","The number of basis functions that is used in the Reduced Basis solve during the Online Stage.");
  params.addRequiredParam<Real>("online_mu", "Current values of the differnt layer for which the RB Method is solved.");
  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the input file. Required for the libMesh functions");
  params.addRequiredParam<std::string>("file_name","");

  return params;
}

DwarfElephantRBSystem::DwarfElephantRBSystem(const InputParameters & params):
  GeneralUserObject(params),
  BlockRestrictable(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _skip_matrix_assembly(false),
  _skip_vector_assembly(false),
  _offline_stage(getParam<bool>("offline_stage")),
  _online_stage(getParam<bool>("online_stage")),
  _store_basis_functions(getParam<bool>("store_basis_functions")),
  _online_N(getParam<unsigned int>("online_N")),
  _online_mu(getParam<Real>("online_mu")),
  _rhs_vector(libmesh_nullptr),
  _residual_non_time_vector(libmesh_nullptr),
  _system_matrix(libmesh_nullptr),
  _system_name(getParam<std::string>("system")),
  _parameters_filename(getParam<std::string>("parameters_filename")),
  _es(_use_displaced ?  _fe_problem.getDisplacedProblem()->es() :  _fe_problem.es()),
  _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _mesh_ptr(&_fe_problem.mesh()),
  _rb_con(_es.add_system<DwarfElephantRBConstruction> ("RBSystem")),
  _file_name(getParam<std::string>("file_name"))
{
}

void
DwarfElephantRBSystem::prepareDataStructuresRB()
{

  // Retrieve the RHS vector for the RB method.
  _rhs_vector = &_sys.get_vector(0);
  _residual_non_time_vector = &_sys.get_vector(1);

  // Retrieve the system matrix for the RB method.
  _system_matrix = &_sys.get_matrix("System Matrix");

//  std::string _system_info = _sys.get_info();
//  _console << _system_info << std::endl;

//  std::string _vec0 = _sys.vector_name(0);
//  std::string _vec1 = _sys.vector_name(1);
//  std::string _vec2 = _sys.vector_name(2);
//  std::string _vec3 = _sys.vector_name(3);
//  std::string _vec4 = _sys.vector_name(4);
//  std::string _vec5 = _sys.vector_name(5);
//  std::string _vec6 = _sys.vector_name(6);
//  std::string _vec7 = _sys.vector_name(7);
//  std::string _vec8 = _sys.vector_name(8);
//
//  _console << _vec0 << std::endl;
//  _console << _vec1 << std::endl;
//  _console << _vec2 << std::endl;
//  _console << _vec3 << std::endl;
//  _console << _vec4 << std::endl;
//  _console << _vec5 << std::endl;
//  _console << _vec6 << std::endl;
//  _console << _vec7 << std::endl;
//  _console << _vec8 << std::endl;

//  /// Print the NumericVectors to file
//  const std::string _suffix = ".txt";
//
//  std::string _file_name_rhs = "RHS"+ _system_name +_suffix;
//  std::string _file_name_residual = "residualNonTime"+ _system_name +_suffix;
//
//  std::filebuf _fb_rhs;
//  std::filebuf _fb_residual;
//
//  _fb_rhs.open(_file_name_rhs,std::ios::out);
//  _fb_residual.open(_file_name_residual,std::ios::out);
//
//  std::ostream _rhs_out (&_fb_rhs);
//  std::ostream _residual_out (&_fb_residual);
//
//  _rhs.print(_rhs_out);
//  _residual_non_time.print(_residual_out);
//
//  /// Print the SparseMatrix to file
//  const std::string _suffix = ".txt";
//
//  std::string _file_name_jacobian = _file_name +_suffix;
//
//  std::filebuf _fb_jacobian;
//
//  _fb_jacobian.open(_file_name_jacobian,std::ios::out);
//
//  std::ostream _jacobian_out (&_fb_jacobian);
//
//  _system_matrix->print(_jacobian_out);
}
void DwarfElephantRBSystem::offlineStage()
{
  // This method performs the offline stage of the RB problem.

  // Get and process the necessary input parameters for the
  // offline stage
  _rb_con.process_parameters_file(_parameters_filename);

  // Print the system informations for the RBConstruction system.
  _rb_con.print_info();

  // Initialize the RB construction. Note, we skip the matrix and vector
  // assembly, since this is already done by MOOSE.
  _rb_con.initialize_rb_construction(/*_skip_matrix_assembly, _skip_vector_assembly*/);


//  std::string _vec0 = _sys.vector_name(0);
//  NumericVector <Number> & _rhs = _rb_con.get_vector(_vec0);
//
//  unsigned int _size_rhs = _rhs.size();
//
//  for (unsigned int i=0; i!=_size_rhs/2; i++)
//  {
//    Number _value = _rhs_vector->el(i);
//    _rhs.set(i,_value);
//    _rhs.set(i+(_size_rhs/2),_value);
//  }

  // Computation of the reduced basis space.
  _rb_con.train_reduced_basis();
//
  // Wrtite the offline data to file (xdr format).
  _rb_con.get_rb_evaluation().legacy_write_offline_data_to_files();

  // If desired, store the basis functions (xdr format).
  if (_store_basis_functions)
  {
    _rb_con.get_rb_evaluation().write_out_basis_functions(_rb_con);
  }

  //
  _rb_con.print_basis_function_orthogonality();


}

void DwarfElephantRBSystem::onlineStage()
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

  // Prepare data import from the Input file for the construction
  // of the RB system.
  GetPot infile (_parameters_filename);

  // Add a new system of the type RBConstruction to the EquationSystems
  // where the RB method is performed.
//  DwarfElephantRBConstruction & _rb_con =
//    _es.add_system<DwarfElephantRBConstruction> ("RBSystem");

  // Intialization of the added system.
  _rb_con.init();

  // Build the RBEvaluation object.
  // Required for both the Offline and Online stage.
  DwarfElephantRBEvaluation _rb_eval(_mesh_ptr->comm());

  // Pass a pointer of the RBEvaluation object to the
  // RBConstruction object
  _rb_con.set_rb_evaluation(_rb_eval);


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
}

void
DwarfElephantRBSystem::execute()
{
  prepareDataStructuresRB();
}

void
DwarfElephantRBSystem::finalize()
{
//  _console << std::endl;
  performRBSystem();
}
