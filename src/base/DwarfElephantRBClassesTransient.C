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

  // Real
  // DwarfElephantRBConstructionTransient::truth_solve(int write_interval)
  // {
  //   LOG_SCOPE("truth_solve()", "TransientRBConstruction");
  //
  //   const RBParameters & mu = get_parameters();
  //   const unsigned int n_time_steps = get_n_time_steps();
  //
  //   //   // NumericVector for computing true L2 error
  //   //   std::unique_ptr<NumericVector<Number>> temp = NumericVector<Number>::build();
  //   //   temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
  //
  //   // Apply initial condition again.
  //   initialize_truth();
  //   set_time_step(0);
  //
  //   // Now compute the truth outputs
  //   for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
  //     {
  //       truth_outputs_all_k[n][0] = 0.;
  //       for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
  //         {
  //           truth_outputs_all_k[n][0] += get_rb_theta_expansion().eval_output_theta(n,q_l,mu)*
  //             get_output_vector(n,q_l)->dot(*solution);
  //         }
  //     }
  //
  //   // Load initial projection error into temporal_data dense matrix
  //   if (compute_truth_projection_error)
  //     set_error_temporal_data();
  //
  //   for (unsigned int time_level=1; time_level<=n_time_steps; time_level++)
  //     {
  //       set_time_step(time_level);
  //
  //       *old_local_solution = *current_local_solution;
  //
  //       // We assume that the truth assembly has been attached to the system
  //       truth_assembly();
  //
  //       // truth_assembly assembles into matrix and rhs, so use those for the solve
  //       solve_for_matrix_and_rhs(*get_linear_solver(), *matrix, *rhs);
  //
  //       // The matrix doesn't change at each timestep, so we
  //       // can set reuse_preconditioner == true
  //       linear_solver->reuse_preconditioner(true);
  //
  //       if (assert_convergence)
  //         {
  //           check_convergence(*get_linear_solver());
  //         }
  //
  //       // Now compute the truth outputs
  //       for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
  //         {
  //           truth_outputs_all_k[n][time_level] = 0.;
  //           for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
  //             {
  //               truth_outputs_all_k[n][time_level] +=
  //                 get_rb_theta_expansion().eval_output_theta(n,q_l,mu)*get_output_vector(n,q_l)->dot(*solution);
  //             }
  //         }
  //
  //       // load projection error into column _k of temporal_data matrix
  //       if (compute_truth_projection_error)
  //         set_error_temporal_data();
  //
  //       if ((write_interval > 0) && (time_level%write_interval == 0))
  //         {
  //           libMesh::out << std::endl << "Truth solve, plotting time step " << time_level << std::endl;
  //
  //           std::ostringstream file_name;
  //
  //           file_name << "truth.e.";
  //           file_name << std::setw(3)
  //                     << std::setprecision(0)
  //                     << std::setfill('0')
  //                     << std::right
  //                     << time_level;
  //
  // #ifdef LIBMESH_HAVE_EXODUS_API
  //           ExodusII_IO(get_mesh()).write_equation_systems (file_name.str(),
  //                                                           this->get_equation_systems());
  // #endif
  //         }
  //     }
  //
  //   // Set reuse_preconditioner back to false for subsequent solves.
  //   linear_solver->reuse_preconditioner(false);
  //
  //   // Get the L2 norm of the truth solution at time-level _K
  //   // Useful for normalizing our true error data
  //   L2_matrix->vector_mult(*inner_product_storage_vector, *solution);
  //   Real final_truth_L2_norm = libmesh_real(std::sqrt(inner_product_storage_vector->dot(*solution)));
  //
  //
  //   return final_truth_L2_norm;
  // }

  // void
  // DwarfElephantRBConstructionTransient::print_info()
  // {
  //   RBConstruction::print_info();
  //
  //   libMesh::out << std::endl << "TransientRBConstruction parameters:" << std::endl;
  //
  //   if (is_rb_eval_initialized())
  //     {
  //       // Print out info that describes the current setup
  //       TransientRBThetaExpansion & trans_theta_expansion =
  //         cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  //       libMesh::out << "Q_m: " << trans_theta_expansion.get_n_M_terms() << std::endl;
  //     }
  //   else
  //     {
  //       libMesh::out << "RBThetaExpansion member is not set yet" << std::endl;
  //     }
  //   libMesh::out << "Number of time-steps: " << get_n_time_steps() << std::endl;
  //   libMesh::out << "dt: " << get_delta_t() << std::endl;
  //   libMesh::out << "euler_theta (time discretization parameter): " << get_euler_theta() << std::endl;
  //   if (get_POD_tol() > 0.)
  //     libMesh::out << "POD_tol: " << get_POD_tol() << std::endl;
  //   if (max_truth_solves > 0)
  //     libMesh::out << "Maximum number of truth solves: " << max_truth_solves << std::endl;
  //   libMesh::out << "delta_N (number of basis functions to add each POD-Greedy step): " << get_delta_N() << std::endl;
  //   if (nonzero_initialization)
  //     {
  //       libMesh::out << "Using initial conditions provided by MOOSE." << std::endl;
  //     }
  //   else
  //     {
  //       libMesh::out << "Using zero initial condition" << std::endl;
  //     }
  //   libMesh::out << std::endl;
  // }
  //
  // void
  // DwarfElephantRBConstructionTransient::initialize_truth()
  // {
  //   if (nonzero_initialization)
  //     {
  //       // DwarfElephantRBEvaluationTransient & trans_rb_eval = cast_ref<DwarfElephantRBEvaluationTransient &>(get_rb_evaluation());
  //       // *this->solution = *trans_rb_eval.get_fe_problem().es().get_system("rb0").solution;
  //       Xdr IC_data(init_filename, READ);
  //       read_serialized_data(IC_data, false);
  //     }
  //   else
  //     {
  //       // Otherwise zero out the solution as a default
  //       this->solution->zero();
  //     }
  //   this->solution->close();
  //   this->update();
  //
  //
  // }

  // Real
  // DwarfElephantRBConstructionTransient::train_reduced_basis(const bool resize_rb_eval_data)
  // {
  //   compute_truth_projection_error = true;
  //   libMesh::out << "Normalize? " << get_normalize_rb_bound_in_greedy () << std::endl;
  //   Real value = train_reduced_basis_steady(resize_rb_eval_data);
  //   compute_truth_projection_error = false;
  //
  //   return value;
  // }

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

  // void
  // DwarfElephantRBConstructionTransient::add_IC_to_RB_space()
  // {
  //   LOG_SCOPE("add_IC_to_RB_space()", "TransientRBConstruction");
  //
  //   if (get_rb_evaluation().get_n_basis_functions() > 0)
  //     libmesh_error_msg("Error: Should not call TransientRBConstruction::add_IC_to_RB_space() " \
  //                       << "on a system that already contains basis functions.");
  //
  //   if (!nonzero_initialization)
  //     libmesh_error_msg("Error: Should not call TransientRBConstruction::add_IC_to_RB_space() " \
  //                       << "when nonzero_initialization==false.");
  //
  //   initialize_truth();
  //
  //   // load the new basis function into the basis_functions vector.
  //   get_rb_evaluation().basis_functions.emplace_back(NumericVector<Number>::build(this->comm()));
  //   NumericVector<Number> & current_bf = get_rb_evaluation().get_basis_function(get_rb_evaluation().get_n_basis_functions()-1);
  //   current_bf.init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
  //   current_bf = *solution;
  //
  //   // We can just set the norm to 1.
  //   get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector,*solution);
  //
  //   Real current_bf_norm = libmesh_real(std::sqrt( current_bf.dot(*inner_product_storage_vector) ));
  //   current_bf.scale(1./current_bf_norm);
  //
  //   unsigned int saved_delta_N = get_delta_N();
  //   set_delta_N(1);
  //   update_system();
  //   set_delta_N(saved_delta_N);
  // }

  // void
  // DwarfElephantRBConstructionTransient::update_system()
  // {
  //   // If delta_N is set to zero, there is nothing to update
  //   if (get_delta_N() == 0)
  //     return;
  //
  //   RBConstruction::update_system();
  //
  //   libMesh::out << "Updating RB initial conditions" << std::endl;
  //   update_RB_initial_condition_all_N();
  // }

  // void
  // DwarfElephantRBConstructionTransient::update_RB_initial_condition_all_N()
  // {
  //   LOG_SCOPE("update_RB_initial_condition_all_N()", "TransientRBConstruction");
  //
  //   libMesh::out << "Start Update RB initial conditions" << std::endl;
  //
  //   TransientRBEvaluation & trans_rb_eval = cast_ref<TransientRBEvaluation &>(get_rb_evaluation());
  //
  //   libMesh::out << "Start Initialize truth" << std::endl;
  //   // Load the initial condition into the solution vector
  //   initialize_truth();
  //   libMesh::out << "End Initialize truth" << std::endl;
  //
  //   std::unique_ptr<NumericVector<Number>> temp1 = NumericVector<Number>::build(this->comm());
  //   temp1->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
  //
  //   std::unique_ptr<NumericVector<Number>> temp2 = NumericVector<Number>::build(this->comm());
  //   temp2->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
  //
  //
  //   unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();
  //
  //   // First compute the right-hand side vector for the L2 projection
  //   L2_matrix->vector_mult(*temp1, *solution);
  //
  //   for (unsigned int i=(RB_size-delta_N); i<RB_size; i++)
  //     {
  //       RB_ic_proj_rhs_all_N(i) = temp1->dot(get_rb_evaluation().get_basis_function(i));
  //     }
  //
  //
  //   // Now compute the projection for each N
  //   DenseMatrix<Number> RB_L2_matrix_N;
  //   DenseVector<Number> RB_rhs_N;
  //   for (unsigned int N=(RB_size-delta_N); N<RB_size; N++)
  //     {
  //       // We have to index here by N+1 since the loop index is zero-based.
  //       trans_rb_eval.RB_L2_matrix.get_principal_submatrix(N+1, RB_L2_matrix_N);
  //
  //       RB_ic_proj_rhs_all_N.get_principal_subvector(N+1, RB_rhs_N);
  //
  //       DenseVector<Number> RB_ic_N(N+1);
  //
  //       // Now solve the linear system
  //       RB_L2_matrix_N.lu_solve(RB_rhs_N, RB_ic_N);
  //
  //       // Load RB_ic_N into RB_initial_condition_all_N
  //       trans_rb_eval.RB_initial_condition_all_N[N] = RB_ic_N;
  //
  //       // Compute the L2 error for the RB initial condition
  //       // This part is dependent on the truth space.
  //
  //       // load the RB solution into temp1
  //       temp1->zero();
  //       for (unsigned int i=0; i<N+1; i++)
  //         {
  //           temp1->add(RB_ic_N(i), get_rb_evaluation().get_basis_function(i));
  //         }
  //
  //       // subtract truth initial condition from RB_ic_N
  //       temp1->add(-1., *solution);
  //
  //       // Compute L2 norm error, i.e. sqrt(M(solution,solution))
  //       temp2->zero();
  //       L2_matrix->vector_mult(*temp2, *temp1);
  //
  //       trans_rb_eval.initial_L2_error_all_N[N] = libmesh_real(std::sqrt(temp2->dot(*temp1)));
  //     }
  // }



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
