///---------------------------------INCLUDE---------------------------------
// Moose includes
#include "FEProblem.h"
#include "FileMesh.h"
#include "MooseApp.h"
#include "NonlinearSystem.h"

// Moose includes (own)
#include "RBOutput.h"

// libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"

#include "libmesh/rb_data_serialization.h"
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
  params.addParam<unsigned int>("test",20,"");

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
    _non_sys_ptr(&_problem_ptr->getNonlinearSystem()),
    _aux_sys_ptr(&_problem_ptr->getAuxiliarySystem()),
    _kernel_warehouse(_non_sys_ptr->getKernelWarehouse()),
    _no_time_kernels_ref(_non_sys_ptr->getNonTimeKernelWarehouse()),
    _mesh_ptr(&_problem_ptr->mesh()),
    _tid(parameters.isParamValid("_tid") ? parameters.get<THREAD_ID>("_tid") : 0),
    _assembly_ptr(&_problem_ptr->assembly(_tid)),
    _test(getParam<unsigned int>("test")),
    _ptr(_problem_ptr->couplingMatrix())
{
}

///-------------------------------------------------------------------------

void
RBOutput::initRBSystem()
{

   GetPot infile (parameters_filename);

  // Create the new EquationSystems
  _es_ptr = new EquationSystems(_mesh_ptr->getMesh());

  SimpleRBConstructionTest & _rb_con =
    _es_ptr->add_system<SimpleRBConstructionTest> ("RBSystem");

  // Intialization of the equation system
  _es_ptr->init();

// Build the RBEvaluation object
// Required for both the Offline and Online stage.
  SimpleRBEvaluation _rb_eval(_mesh_ptr->comm());
//
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

//    Computation of the reduced basis space by taking snapshots.
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

//void
//test::interior_assembly()
//{}

void
RBOutput::output(const ExecFlagType & type)
{
//KernelBase * k = dynamic_cast<KernelBase *>(this);
  if (type==EXEC_FINAL)
  {

    _console << std::endl;
    _test_sys = _non_sys_ptr->_num_residual_evaluations;
    MooseSharedPointer<MooseObject> _moose_object_ptr;
//    _block_test = _jacobian_block_ptr->_ivar;
//    _jacobian_test = &_jacobian_block_ptr->_jacobian;
//    const MooseArray<Point> & _test_test = _assembly_ptr->normals();
//    _jacobian = _rb_ptr->_local_ke;

    //_assembly_ptr->cacheJacobian();
//    _console << _mesh_ptr << std::endl;
//    _console << _non_sys_ptr << std::endl;
//    _console << _assembly_ptr << std::endl;
//     initRBSystem();
  }
    if (type==EXEC_TIMESTEP_END)
  {
    FEProblem & _problem = *_problem_ptr;
    const std::string _name_test = "RB";
    const std::string & _name_test_ref = _name_test;
    test(_problem,_name_test_ref);
    _console << std::endl;
    _console << _aux_sys_ptr << std::endl;
    _console << _name_test << std::endl;
    _console << _name_test_ref << std::endl;
    _console << std::endl;
    initRBSystem();


  }
}
