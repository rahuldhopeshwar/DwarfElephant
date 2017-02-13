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
  params.addRequiredParam<bool>("offline_stage", "Determines whether the \
                                Offline stage will be calculated or not.");
  params.addRequiredParam<bool>("online_stage", "Determines whether the \
                                Online stage will be calculated or not.");
  params.addRequiredParam<bool>("store_basis_functions","Determines whether \
                                the basis functions are stored or not");
  params.addRequiredParam<unsigned int>("online_N","The number of basis \
                                        functions that is used in the \
                                        Reduced Basis solve during the \
                                        Online Stage.");
  params.addRequiredParam<Real>("online_mu0", "Current value for which the \
                                 RB Method is solved.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the \
                                        input file. Required for the libMesh\
                                        functions");
//  params.addRequiredParam<NonlinearVariableName>("variable", "The name of \
                                                  the variable that this \
                                                  class should output.");

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
    _online_mu0_parameters(getParam<Real>("online_mu0")),
    _mesh_ptr(&_problem_ptr->mesh())
//    _tid(parameters.get<THREAD_ID>("_tid")),
////    _non_sys_ptr(&_problem_ptr->getNonlinearSystem()),
//    _aux_sys_ptr(&_problem_ptr->getAuxiliarySystem()),
//    _n_aux_var(_aux_sys_ptr->nVariables()),
////    _var(_aux_sys_ptr->getVariable(_tid, parameters.get<NonlinearVariableName>("variable")))
//    _var(_aux_sys_ptr->getVariable(_tid, 1))
//    _kernel_warehouse(_non_sys_ptr->getKernelWarehouse()),
//    _no_time_kernels_ref(_non_sys_ptr->getNonTimeKernelWarehouse()),
//    _assembly_ptr(&_problem_ptr->assembly(_tid)),
//    _test(getParam<unsigned int>("test")),
//    _ptr(_problem_ptr->couplingMatrix()),
//    _non_solver_ptr(_non_sys_ptr->nonlinearSolver()),
//    _phi_test(_assembly_ptr->phi()),
//    _executioner(_app.executioner())
//    kernels(_kernel_warehouse.getObjects(_tid))

{
}

///-------------------------------------------------------------------------

void
RBOutput::performRBSystem()
{

  GetPot infile (parameters_filename);

  // Create the new EquationSystems
  _es_ptr = new EquationSystems(_mesh_ptr->getMesh());

  RBSimpleConstruction & _rb_con =
    _es_ptr->add_system<RBSimpleConstruction> ("RBSystem");

  // Intialization of the equation system
  _es_ptr->init();

  // Build the RBEvaluation object
  // Required for both the Offline and Online stage.
  RBSimpleEvaluation _rb_eval(_mesh_ptr->comm());

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

    // Computation of the reduced basis space by taking snapshots.
    _rb_con.train_reduced_basis();


    _rb_con.get_rb_evaluation().legacy_write_offline_data_to_files();

    if (_store_basis_functions)
    {
      _rb_con.get_rb_evaluation().write_out_basis_functions(_rb_con);
    }

    _rb_con.print_basis_function_orthogonality();

//    /// Online stage
//    RBParameters _online_mu;
//
//    _online_mu.set_value("mu0", _online_mu0_parameters);
//    _rb_eval.set_parameters(_online_mu);
//    _rb_eval.rb_solve(_online_N);
//
//    _rb_eval.print_parameters();
   }
//////////    elseif(_offline_stage==false && _online_stage==true)
//////////     _rb_eval.legacy_read_offline_data_from_files();
//////////    elseif(_offline_stage==true && _online_stage==false)
//////////    elseif(_offline_stage==false && _online_stage==false
////////    //return _rb_con;
  }

void
RBOutput::output(const ExecFlagType & type)
{
  if (type==EXEC_TIMESTEP_END)
  {
//    std::ofstream _file_output;
//    _file_output.open("test3.txt");
//
//    const NumericVector< Number > * _residual_current = _aux_sys_ptr->currentSolution();
//    NumericVector< Number > & _residual = _aux_sys_ptr->solution();
//    NumericVector< Number > & _residual_old = _aux_sys_ptr->solutionOld();
//    NumericVector< Number > & _residual_older = _aux_sys_ptr->solutionOlder();
//
//    _file_output << _residual_current << std::endl;
////    _file_output << std::endl;
//    _file_output << _residual << std::endl;
////    _file_output << std::endl;
//    _file_output << _residual_old << std::endl;
////    _file_output << std::endl;
//    _file_output << _residual_older << std::endl;
//    _file_output.close();
//
//    _console << std::endl;
//    _console << _n_aux_var << std::endl;
//    _console << std::endl;
//    _console << _residual << std::endl;
//    _console << std::endl;

//    performRBSystem();
  }
}
