/**
 * In this class simplified subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
 *
 * RBSimpleEvaluation: requires only the definition of the lower coercivity
 * constant. The value is here specified for a Conduction problem.
 *
 * RBSimpleConstruction: In order to construct the RB System with the
 * RBSimpleEvaluation subclass the method build_rb_evaluation needs to be
 * overriden.
 *
 * NOTE: ENSURE THAT THE CLASS IS INHERITING FROM THE CORRECT RBSTRUCTURES.
 */

///-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
#include "DwarfElephantRBClassesTransient.h"

DwarfElephantRBConstructionTransient::DwarfElephantRBConstructionTransient (EquationSystems & es,
                      const std::string & name_in,
                      const unsigned int number_in)
  : Parent(es, name_in, number_in)
{}

void
DwarfElephantRBConstructionTransient::init_data()
  {
    u_var = this->add_variable (get_equation_systems().get_system(0).variable_name(0) + "(RB)", libMesh::FIRST);

    Parent::init_data();
  }

  void
  DwarfElephantRBConstructionTransient::print_info()
  {
    RBConstruction::print_info();

    libMesh::out << std::endl << "TransientRBConstruction parameters:" << std::endl;

    if (is_rb_eval_initialized())
      {
        // Print out info that describes the current setup
        TransientRBThetaExpansion & trans_theta_expansion =
          cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
        libMesh::out << "Q_m: " << trans_theta_expansion.get_n_M_terms() << std::endl;
      }
    else
      {
        libMesh::out << "RBThetaExpansion member is not set yet" << std::endl;
      }
    libMesh::out << "Number of time-steps: " << get_n_time_steps() << std::endl;
    libMesh::out << "dt: " << get_delta_t() << std::endl;
    libMesh::out << "euler_theta (time discretization parameter): " << get_euler_theta() << std::endl;
    if (get_POD_tol() > 0.)
      libMesh::out << "POD_tol: " << get_POD_tol() << std::endl;
    if (max_truth_solves > 0)
      libMesh::out << "Maximum number of truth solves: " << max_truth_solves << std::endl;
    libMesh::out << "delta_N (number of basis functions to add each POD-Greedy step): " << get_delta_N() << std::endl;
    if (nonzero_initialization)
      {
        libMesh::out << "Using initial conditions provided by MOOSE." << std::endl;
      }
    else
      {
        libMesh::out << "Using zero initial condition" << std::endl;
      }
    libMesh::out << std::endl;
  }

  void
  DwarfElephantRBConstructionTransient::initialize_truth()
  {
    if (nonzero_initialization)
      {
        DwarfElephantRBEvaluationTransient & trans_rb_eval = cast_ref<DwarfElephantRBEvaluationTransient &>(get_rb_evaluation());
        *this->solution.get() = *trans_rb_eval.get_fe_problem().es().get_system("rb0").solution.get();
      }
    else
      {
        // Otherwise zero out the solution as a default
        this->solution->zero();
      }
    this->solution->close();
    this->update();
  }

  Real DwarfElephantRBConstructionTransient::get_RB_error_bound()
  {
    get_rb_evaluation().set_parameters( get_parameters() );

    Real error_bound = get_rb_evaluation().rb_solve(get_rb_evaluation().get_n_basis_functions());

    if (get_normalize_rb_bound_in_greedy())
      {
        Real error_bound_normalization = get_rb_evaluation().rb_solve(0);

        if ((error_bound < get_abs_training_tolerance ()) ||
            (error_bound_normalization < get_abs_training_tolerance ()))
          {
            // We don't want to normalize this error bound if the bound or the
            // normalization value are below the absolute tolerance. Hence do nothing
            // in this case.
          }
        else
          error_bound /= error_bound_normalization;
      }

    return error_bound;
  }


///------------------------DWARFELEPHANTRBEVALUATION------------------------
  DwarfElephantRBEvaluationTransient::DwarfElephantRBEvaluationTransient(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem):
    TransientRBEvaluation(comm),
    fe_problem(fe_problem)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }

  Real
  DwarfElephantRBEvaluationTransient::get_stability_lower_bound()
  {
    const RBParameters & mu = get_parameters();
    bool norm_values = fe_problem.getUserObject<DwarfElephantOfflineOnlineStageTransient>("performRBSystem")._norm_online_values;
    unsigned int norm_id = fe_problem.getUserObject<DwarfElephantOfflineOnlineStageTransient>("performRBSystem")._norm_id;

    Real min_mu;
    Real min_mu_i;

    min_mu = mu.get_value("mu_0");

    if(norm_values)
      min_mu = min_mu/mu.get_value("mu_"+ std::to_string(norm_id));

    for (unsigned int  i = 1; i != mu.n_parameters(); i++)
    {
      const std::string mu_name = "mu_" + std::to_string(i);
      if(norm_values)
        min_mu_i = std::min(min_mu, mu.get_value(mu_name)/mu.get_value(mu_name));
      else
        min_mu_i = std::min(min_mu, mu.get_value(mu_name));

      if (min_mu_i < min_mu)
        min_mu = min_mu_i;
    }

    return min_mu;
  }
