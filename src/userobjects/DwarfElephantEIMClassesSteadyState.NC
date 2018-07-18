#include "DwarfElephantEIMClassesSteadyState.h"

// A simple subclass of RBEIMEvaluation. Overload
// evaluate_parametrized_function to define the
// function that we "empirically" interpolate.
DwarfElephantEIMEvaluationSteadyState::DwarfElephantEIMEvaluationSteadyState(const libMesh::Parallel::Communicator & comm) :
    RBEIMEvaluation(comm)
  {
    attach_parametrized_function(&sg);
  }


  /**
   * Constructor.
   */
  DwarfElephantEIMConstructionSteadyState::DwarfElephantEIMConstructionSteadyState (EquationSystems & es,
                         const std::string & name_in,
                         const unsigned int number_in)
    : Parent(es, name_in, number_in)
  {
  }

  /**
   * Initialize the implicit system that is used to perform L2 projections.
   */
  void DwarfElephantEIMConstructionSteadyState::init_implicit_system()
  {
    this->add_variable ("L2_proj_var", libMesh::FIRST);
  }

  /**
   * Initialize the explicit system that is used to store the basis functions.
   */
  void DwarfElephantEIMConstructionSteadyState::init_explicit_system()
  {
    u_var = get_explicit_system().add_variable ("f_EIM", libMesh::FIRST);
  }

