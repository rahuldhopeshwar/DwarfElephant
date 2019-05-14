#include "DwarfElephantRBIterationAdaptiveDT.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBIterationAdaptiveDT);

template <>
InputParameters
validParams<DwarfElephantRBIterationAdaptiveDT>()
{
  InputParameters params = validParams<DwarfElephantRBExecutioner>();
  params.addClassDescription("Adjust the timestep based on the number of iterations");
  params.addParam<int>("optimal_iterations",
                       "The target number of nonlinear iterations for adaptive timestepping");
  params.addParam<int>("iteration_window",
                       "Attempt to grow/shrink timestep if the iteration count "
                       "is below/above 'optimal_iterations plus/minus "
                       "iteration_window' (default = optimal_iterations/5).");
  params.addParam<unsigned>("linear_iteration_ratio",
                            "The ratio of linear to nonlinear iterations "
                            "to determine target linear iterations and "
                            "window for adaptive timestepping (default = "
                            "25)");
  params.addParam<Real>("growth_factor",
                        2.0,
                        "Factor to apply to timestep if easy convergence (if "
                        "'optimal_iterations' is specified) or if recovering "
                        "from failed solve");
  params.addParam<Real>("cutback_factor",
                        0.5,
                        "Factor to apply to timestep if difficult "
                        "convergence (if 'optimal_iterations' is specified) "
                        "or if solution failed");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");

  params.declareControllable("growth_factor cutback_factor");
  return params;
}

DwarfElephantRBIterationAdaptiveDT::DwarfElephantRBIterationAdaptiveDT(const InputParameters & parameters) :
 DwarfElephantRBExecutioner(parameters),
 _linear_iteration_ratio(isParamValid("linear_iteration_ratio")
                             ? getParam<unsigned>("linear_iteration_ratio")
                             : 25), // Default to 25
 _growth_factor(getParam<Real>("growth_factor")),
 _cutback_factor(getParam<Real>("cutback_factor")),
 _cutback_occurred(declareRestartableData<bool>("cutback_occurred", false))
{}

void
DwarfElephantRBIterationAdaptiveDT::execute()
{
  DwarfElephantRBExecutioner::execute();
}

Real
DwarfElephantRBIterationAdaptiveDT::computeDT()
{
  UserObjectName _initialize_rb_system_name = getParam<UserObjectName>("initial_rb_userobject");
  const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = _fe_problem.getUserObject<DwarfElephantInitializeRBSystemTransient>(_initialize_rb_system_name);
  //
  const unsigned int growth_l_its(_optimal_iterations > _iteration_window
                                      ? _linear_iteration_ratio *
                                      (_optimal_iterations - _iteration_window)
                                      : 0);
  const unsigned int shrink_l_its(_linear_iteration_ratio *
                                    (_optimal_iterations + _iteration_window));

  dt = _initialize_rb_system._rb_con_ptr->get_delta_t();
  _n_l_iter = _initialize_rb_system._rb_con_ptr->n_linear_iterations();

  // Adaptive Time stepping
  // TODO: Add further functionality
  if(_n_l_iter < growth_l_its)
    dt *= _growth_factor;
  else if (_n_l_iter > shrink_l_its)
    dt *= _cutback_factor;

  return dt;
}
//
