/**
 * In this class subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
 *
 * NOTE: ENSURE THAT THE CLASS IS USING THE CORRECT RBSTRUCTURES.
 */

 //-------------------------------------------------------------------------
 #include "DwarfElephantRBClassesSteadyState.h"

DwarfElephantRBConstructionSteadyState::DwarfElephantRBConstructionSteadyState (EquationSystems & es,
                                                                                const std::string & name_in,
                                                                                const unsigned int number_in)
    : Parent(es, name_in, number_in)
  {}

void
DwarfElephantRBConstructionSteadyState::init_data()
{

  const std::set<subdomain_id_type> lm_subdomains {11};
  u_var = this->add_variable(get_equation_systems().get_system(0).variable_name(0) + "(RB)");
  lm_var = this->add_variable(get_equation_systems().get_system(0).variable_name(1) + "(RB)", FIRST, LAGRANGE, &lm_subdomains);

  //get_equation_systems().update();
  Parent::init_data();
  get_equation_systems().print_info();
}

Real
DwarfElephantRBConstructionSteadyState::truth_solve(int plot_solution)
{
LOG_SCOPE("truth_solve()", "RBConstruction");

truth_assembly();

// truth_assembly assembles into matrix and rhs, so use those for the solve
if (extra_linear_solver)
  {
    // If extra_linear_solver has been initialized, then we use it for the
    // truth solves.
    solve_for_matrix_and_rhs(*extra_linear_solver, *matrix, *rhs);

    if (assert_convergence)
      check_convergence(*extra_linear_solver);
  }
else
  {
    solve_for_matrix_and_rhs(*get_linear_solver(), *matrix, *rhs);

    if (assert_convergence)
      check_convergence(*get_linear_solver());
  }



const RBParameters & mu = get_parameters();

for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
  {
    truth_outputs[n] = 0.;
    for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
      truth_outputs[n] += get_rb_theta_expansion().eval_output_theta(n, q_l, mu)*
        get_output_vector(n,q_l)->dot(*solution);
  }

if (plot_solution > 0)
  {
    // or set solution ?
    // this->get_equation_systems().solution = solution;
#if defined(LIBMESH_USE_COMPLEX_NUMBERS)
    GMVIO(get_mesh()).write_equation_systems ("truth.gmv",
                                              this->get_equation_systems());
#else
#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(get_mesh()).write_equation_systems ("truth.e",
                                                        this->get_equation_systems());
#endif
#endif

//#if defined(LIBMESH_HAVE_VTK)
  //VTK_IO(get_mesh()).write ("mesh_truth.pvtu");
  //VTK_IO(get_mesh()).write_equation_systems ("truth.pvtu", this->get_equation_systems());
//#endif


  }

  get_vector(0).zero();
  get_vector(0).add(*solution->clone());

ExodusII_IO(get_mesh()).write_equation_systems ("truth.e",
                                                    this->get_equation_systems());

                                                 std::filebuf _buffer;
                                                  _buffer.open("solution_dump.dat", std::ios::app);
                                                    std::ostream state_file(&_buffer);
                                                    solution->print(state_file);
                                        _buffer.close();
// Get the X norm of the truth solution
// Useful for normalizing our true error data
inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
Number truth_X_norm = std::sqrt(inner_product_storage_vector->dot(*solution));
std::cout<< "solution norm "<< truth_X_norm << std::endl;
return libmesh_real(truth_X_norm);
}

Real
DwarfElephantRBConstructionSteadyState::custom_train_reduced_basis(const bool resize_rb_eval_data)
{
  libMesh::out << "custom training = "
              << std::endl;

 LOG_SCOPE("train_reduced_basis()", "RBConstruction");

 int count = 0;

 // initialize rb_eval's parameters
 get_rb_evaluation().initialize_parameters(*this);

 // possibly resize data structures according to Nmax
 if (resize_rb_eval_data)
   {
     get_rb_evaluation().resize_data_structures(get_Nmax());
   }

 // Clear the Greedy param list
 for (std::size_t i=0; i<get_rb_evaluation().greedy_param_list.size(); i++)
   get_rb_evaluation().greedy_param_list[i].clear();

 get_rb_evaluation().greedy_param_list.clear();

 Real training_greedy_error;


 // If we are continuing from a previous training run,
 // we might already be at the max number of basis functions.
 // If so, we can just return.
 if (get_rb_evaluation().get_n_basis_functions() >= get_Nmax())
   {
     libMesh::out << "Maximum number of basis functions reached: Nmax = "
                  << get_Nmax() << std::endl;
     return 0.;
   }


 // Compute the dual norms of the outputs if we haven't already done so
 compute_output_dual_innerprods();

inner_product_matrix->print_matlab();
 // Compute the Fq Riesz representor dual norms if we haven't already done so
 compute_Fq_representor_innerprods();

 libMesh::out << std::endl << "---- Performing Greedy basis enrichment ----" << std::endl;
 Real initial_greedy_error = 0.;
 bool initial_greedy_error_initialized = false;
 while (true)
   {
     libMesh::out << std::endl << "---- Basis dimension: "
                  << get_rb_evaluation().get_n_basis_functions() << " ----" << std::endl;

     if (count > 0 || (count==0 && use_empty_rb_solve_in_greedy))
       {
         libMesh::out << "Performing RB solves on training set" << std::endl;
         training_greedy_error = compute_max_error_bound();

         libMesh::out << "Maximum error bound is " << training_greedy_error << std::endl << std::endl;

         // record the initial error
         if (!initial_greedy_error_initialized)
           {
             initial_greedy_error = training_greedy_error;
             initial_greedy_error_initialized = true;
           }

         // Break out of training phase if we have reached Nmax
         // or if the training_tolerance is satisfied.
         if (greedy_termination_test(training_greedy_error, initial_greedy_error, count))
           break;
       }

     libMesh::out << "Performing truth solve at parameter:" << std::endl;
     print_parameters();

     // Update the list of Greedily selected parameters
     this->update_greedy_param_list();

     // Perform an Offline truth solve for the current parameter
     Real truth_X_norm = truth_solve(1);
     libMesh::out << truth_X_norm << std::endl;

     // Add orthogonal part of the snapshot to the RB space
     libMesh::out << "Enriching the RB space" << std::endl;
     enrich_RB_space();

     update_system();

     // Increment counter
     count++;
   }
 this->update_greedy_param_list();

 return training_greedy_error;
}

 Real
 DwarfElephantRBConstructionSteadyState::compute_residual_dual_norm(const unsigned int N)
{
   LOG_SCOPE("compute_residual_dual_norm()", "RBConstruction");

   // Put the residual in rhs in order to compute the norm of the Riesz representor
   // Note that this only works in serial since otherwise each processor will
   // have a different parameter value during the Greedy training.

   UniquePtr< NumericVector<Number> > RB_sol = NumericVector<Number>::build(this->comm());
   RB_sol->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

   UniquePtr< NumericVector<Number> > temp = NumericVector<Number>::build(this->comm());
   temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

   for(unsigned int i=0; i<N; i++)
   {
     RB_sol->add(get_rb_evaluation().RB_solution(i), get_rb_evaluation().get_basis_function(i));
   }

   this->truth_assembly();
   matrix->vector_mult(*temp, *RB_sol);
   rhs->add(-1., *temp);

   // Then solve to get the Reisz representor
   matrix->zero();
   matrix->add(1., *inner_product_matrix);
//   if(constrained_problem)
//     matrix->add(1., *constraint_matrix);

   solution->zero();
   solve();
   // Make sure we didn't max out the number of iterations
   if( (this->n_linear_iterations() >=
        this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
       (this->final_linear_residual() >
        this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
   {
     libmesh_error_msg("Warning: Linear solver may not have converged! Final linear residual = "
                       << this->final_linear_residual() << ", number of iterations = "
                       << this->n_linear_iterations());
   }

   inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);

   Real slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);

   return std::sqrt( libmesh_real(slow_residual_norm_sq) );
}

DwarfElephantRBEvaluationSteadyState::DwarfElephantRBEvaluationSteadyState(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem):
    RBEvaluation(comm),
    fe_problem(fe_problem)
{
  set_rb_theta_expansion(_rb_theta_expansion);
}

Real
DwarfElephantRBEvaluationSteadyState::get_stability_lower_bound()
{
  const RBParameters & mu = get_parameters();

  Real min_mu = mu.get_value("mu_0");

  for (unsigned int  i = 1; i != mu.n_parameters(); i++)
  {
    const std::string mu_name = "mu_" + std::to_string(i);
    Real min_mu_i = std::min(min_mu, mu.get_value(mu_name));

    if (min_mu_i < min_mu)
      min_mu = min_mu_i;
  }

  return min_mu;
}
