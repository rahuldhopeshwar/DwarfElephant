 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemAction.h"

template<>
InputParameters validParams<DwarfElephantInitializeRBSystemAction>()
{
  InputParameters params = validParams<Action>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
  params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
  params.addRequiredParam<std::string>("parameters_filename","Path to the input file. Required for the libMesh functions");

  return params;
}

DwarfElephantRBSystem::DwarfElephantRBSystem(InputParameters params):
  Action(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _parameters_filename(getParam<std::string>("parameters_filename")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _mesh_ptr(&_fe_problem.mesh())
{
}

void
DwarfElephantInitializeRBSystemAction::initializeRBSystem()
{
  // Define the parameter file for the libMesh functions.
  GetPot infile (_parameters_filename);

  // Add a new equation system for the RB construction.
  _rb_con_ptr = &_es.add_system<DwarfElephantRBConstruction> ("RBSystem");

  // Initialization of the added equation system
  _rb_con_ptr->init();

  // Build the RBEvaluation object
  // Required for both the Offline and Online stage.
  DwarfElephantRBEvaluation _rb_eval(_mesh_ptr->comm());

  // Pass a pointer of the RBEvaluation object to the
  // RBConstruction object
  _rb_con_ptr->set_rb_evaluation(_rb_eval);

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
DwarfElephantInitializeRBSystemAction::act()
{
    initializeRBSystem();
}
