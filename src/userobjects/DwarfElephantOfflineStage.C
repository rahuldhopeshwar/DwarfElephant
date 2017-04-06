///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantOfflineStage.h"

template<>
InputParameters validParams<DwarfElephantOfflineStage>()
{
    InputParameters params = validParams<GeneralUserObject>();

    params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
    params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not.");
    params.addParam<bool>("compliant", true, "Determines whether F is equal to the output vector or not.");
    params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
    params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
    params.addParam<bool>("online_stage", false, "Determines whether the Online stage will be calculated or not.");
    params.addParam<std::string>("system","nl0","The name of the system that should be read in.");
    params.addRequiredParam<std::string>("residual_name","Name of the residual vector that is retrieved. The name is either '_Re_non_time' or '_Re_time'.");
    params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");
    params.addParam<Real>("mu_bar", 1., "Value for mu-bar");
    params.addRequiredParam<unsigned int>("online_N","The number of basis functions that is used in the Reduced Basis solve during the Online Stage.");
    params.addRequiredParam<std::vector<Real>>("online_mu", "Current values of the different layers for which the RB Method is solved.");
    params.addRequiredParam<FunctionName>("cache_stiffness_matrix", "");

    return params;
}

DwarfElephantOfflineStage::DwarfElephantOfflineStage(const InputParameters & params):
    GeneralUserObject(params),
    _use_displaced(getParam<bool>("use_displaced")),
    _store_basis_functions(getParam<bool>("store_basis_functions")),
    _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
    _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
    _compliant(getParam<bool>("compliant")),
    _online_stage(getParam<bool>("online_stage")),
    _system_name(getParam<std::string>("system")),
    _residual_name(getParam<std::string>("residual_name")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
    _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject")),
    _function(&getFunction("cache_stiffness_matrix")),
    _mesh_ptr(&_fe_problem.mesh()),
    _subdomain_ids(_mesh_ptr->meshSubdomains()),
    _mu_bar(getParam<Real>("mu_bar")),
    _online_N(getParam<unsigned int>("online_N")),
    _online_mu_parameters(getParam<std::vector<Real>>("online_mu")),
    _nodal_bcs(false)
{
  _cache_stiffness_matrix = dynamic_cast<CacheStiffnessMatrix *>(_function);
}

Real
DwarfElephantOfflineStage::trainReducedBasis(const bool resize_rb_eval_data)
{
  LOG_SCOPE("train_reduced_basis()", "RBConstruction");

  int count = 0;

  // initialize rb_eval's parameters
  _initialize_rb_system._rb_con_ptr->get_rb_evaluation().initialize_parameters(*_initialize_rb_system._rb_con_ptr);

  // possibly resize data structures according to Nmax
  if(resize_rb_eval_data)
    {
      _initialize_rb_system._rb_con_ptr->get_rb_evaluation().resize_data_structures(_initialize_rb_system._rb_con_ptr->get_Nmax());
    }

  // Clear the Greedy param list
  for (std::size_t i=0; i<_initialize_rb_system._rb_con_ptr->get_rb_evaluation().greedy_param_list.size(); i++)
    _initialize_rb_system._rb_con_ptr->get_rb_evaluation().greedy_param_list[i].clear();

  _initialize_rb_system._rb_con_ptr->get_rb_evaluation().greedy_param_list.clear();

  Real training_greedy_error;


  // If we are continuing from a previous training run,
  // we might already be at the max number of basis functions.
  // If so, we can just return.
  if(_initialize_rb_system._rb_con_ptr->get_rb_evaluation().get_n_basis_functions() >= _initialize_rb_system._rb_con_ptr->get_Nmax())
    {
      libMesh::out << "Maximum number of basis functions reached: Nmax = "
                   << _initialize_rb_system._rb_con_ptr->get_Nmax() << std::endl;
      return 0.;
    }


  // Compute the dual norms of the outputs if we haven't already done so
  _initialize_rb_system._rb_con_ptr->compute_output_dual_innerprods();

  // Compute the Fq Riesz representor dual norms if we haven't already done so
  _initialize_rb_system._rb_con_ptr->compute_Fq_representor_innerprods();

  libMesh::out << std::endl << "---- Performing Greedy basis enrichment ----" << std::endl;
  Real initial_greedy_error = 0.;
  bool initial_greedy_error_initialized = false;
  while(true)
    {
      libMesh::out << std::endl << "---- Basis dimension: "
                   << _initialize_rb_system._rb_con_ptr->get_rb_evaluation().get_n_basis_functions() << " ----" << std::endl;

      if( count > 0 || (count==0 && _initialize_rb_system._rb_con_ptr->use_empty_rb_solve_in_greedy) )
        {
          libMesh::out << "Performing RB solves on training set" << std::endl;
          training_greedy_error = _initialize_rb_system._rb_con_ptr->compute_max_error_bound();

          libMesh::out << "Maximum error bound is " << training_greedy_error << std::endl << std::endl;

          // record the initial error
          if (!initial_greedy_error_initialized)
            {
              initial_greedy_error = training_greedy_error;
              initial_greedy_error_initialized = true;
            }

          // Break out of training phase if we have reached Nmax
          // or if the training_tolerance is satisfied.
          if (_initialize_rb_system._rb_con_ptr->greedy_termination_test(training_greedy_error, initial_greedy_error, count))
            break;
        }

      libMesh::out << "Performing truth solve at parameter:" << std::endl;
      _initialize_rb_system._rb_con_ptr->print_parameters();

      // Update the list of Greedily selected parameters
      _initialize_rb_system._rb_con_ptr->update_greedy_param_list();

      // Perform an Offline truth solve for the current parameter
      truthSolve(-1);

      // Add orthogonal part of the snapshot to the RB space
      libMesh::out << "Enriching the RB space" << std::endl;
      _initialize_rb_system._rb_con_ptr->enrich_RB_space();

      _initialize_rb_system._rb_con_ptr->update_system();

      // Increment counter
      count++;
    }
  _initialize_rb_system._rb_con_ptr->update_greedy_param_list();

  return training_greedy_error;
}

Real
DwarfElephantOfflineStage::truthSolve(int plot_solution)
{
  LOG_SCOPE("truth_solve()", "RBConstruction");

  truthAssembly();

  // truth_assembly assembles into matrix and rhs, so use those for the solve
  if (_initialize_rb_system._rb_con_ptr->extra_linear_solver)
    {
      // If extra_linear_solver has been initialized, then we use it for the
      // truth solves.
      _initialize_rb_system._rb_con_ptr->solve_for_matrix_and_rhs(*_initialize_rb_system._rb_con_ptr->extra_linear_solver, *_initialize_rb_system._rb_con_ptr->matrix, *_initialize_rb_system._rb_con_ptr->rhs);

      if (_initialize_rb_system._rb_con_ptr->assert_convergence)
        _initialize_rb_system._rb_con_ptr->check_convergence(*_initialize_rb_system._rb_con_ptr->extra_linear_solver);
    }
  else
    {
      _initialize_rb_system._rb_con_ptr->solve_for_matrix_and_rhs(*_initialize_rb_system._rb_con_ptr->get_linear_solver(), *_initialize_rb_system._rb_con_ptr->matrix, *_initialize_rb_system._rb_con_ptr->rhs);

      if (_initialize_rb_system._rb_con_ptr->assert_convergence)
        _initialize_rb_system._rb_con_ptr->check_convergence(*_initialize_rb_system._rb_con_ptr->get_linear_solver());
    }



  const RBParameters & mu = _initialize_rb_system._rb_con_ptr->get_parameters();

  for(unsigned int n=0; n<_initialize_rb_system._rb_con_ptr->get_rb_theta_expansion().get_n_outputs(); n++)
    {
      _initialize_rb_system._rb_con_ptr->truth_outputs[n] = 0.;
      for(unsigned int q_l=0; q_l<_initialize_rb_system._rb_con_ptr->get_rb_theta_expansion().get_n_output_terms(n); q_l++)
        _initialize_rb_system._rb_con_ptr->truth_outputs[n] += _initialize_rb_system._rb_con_ptr->get_rb_theta_expansion().eval_output_theta(n, q_l, mu)*
          _initialize_rb_system._rb_con_ptr->get_output_vector(n,q_l)->dot(*_initialize_rb_system._rb_con_ptr->solution);
    }

  if(plot_solution > 0)
    {
#if defined(LIBMESH_USE_COMPLEX_NUMBERS)
      GMVIO(_initialize_rb_system._rb_con_ptr->get_mesh()).write_equation_systems ("truth.gmv",
                                                _initialize_rb_system._rb_con_ptr->get_equation_systems());
#else
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(_initialize_rb_system._rb_con_ptr->get_mesh()).write_equation_systems ("truth.e",
                                                      _initialize_rb_system._rb_con_ptr->get_equation_systems());
#endif
#endif
    }

  // Get the X norm of the truth solution
  // Useful for normalizing our true error data
  _initialize_rb_system._rb_con_ptr->inner_product_matrix->vector_mult(*_initialize_rb_system._rb_con_ptr->inner_product_storage_vector, *_initialize_rb_system._rb_con_ptr->solution);
  Number truth_X_norm = std::sqrt(_initialize_rb_system._rb_con_ptr->inner_product_storage_vector->dot(*_initialize_rb_system._rb_con_ptr->solution));

  return libmesh_real(truth_X_norm);
}

  void
  DwarfElephantOfflineStage::truthAssembly()
{
  LOG_SCOPE("truth_assembly()", "RBConstruction");

  const RBParameters & mu = _initialize_rb_system._rb_con_ptr->get_parameters();

  _initialize_rb_system._rb_con_ptr->matrix->zero();
  _initialize_rb_system._rb_con_ptr->rhs->zero();

  _initialize_rb_system._rb_con_ptr->matrix->close();
  _initialize_rb_system._rb_con_ptr->rhs->close();

  {
    // We should have already assembled the matrices
    // and vectors in the affine expansion, so
    // just use them

    for(unsigned int q_a=0; q_a<_initialize_rb_system._rb_con_ptr->get_rb_theta_expansion().get_n_A_terms(); q_a++)
      {
          _initialize_rb_system._rb_con_ptr->matrix->add(_initialize_rb_system._rb_con_ptr->get_rb_theta_expansion().eval_A_theta(q_a, mu), *_initialize_rb_system._rb_con_ptr->get_Aq(q_a));
          _console << _initialize_rb_system._rb_con_ptr->get_rb_theta_expansion().eval_A_theta(q_a, mu) << std::endl;
//          _console << *_initialize_rb_system._rb_con_ptr->matrix << std::endl;
//        _initialize_rb_system._rb_con_ptr->matrix->add(_mu.get_value(_mu_name), *_initialize_rb_system._rb_con_ptr->get_Aq(q_a));
//        PetscMatrix<Number> * _petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_initialize_rb_system._rb_con_ptr->matrix);
//        MatSetOption(_petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//
//        _cache_stiffness_matrix->setCachedSubdomainStiffnessMatrixContributions(*_initialize_rb_system._rb_con_ptr->matrix, q_a);
//      _initialize_rb_system._jacobian_subdomain[*it] ->close();

      }

    UniquePtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build(_initialize_rb_system._rb_con_ptr->comm());
    temp_vec->init (_initialize_rb_system._rb_con_ptr->n_dofs(), _initialize_rb_system._rb_con_ptr->n_local_dofs(), false, PARALLEL);
    for(unsigned int q_f=0; q_f<_initialize_rb_system._rb_con_ptr->get_rb_theta_expansion().get_n_F_terms(); q_f++)
      {
        *temp_vec = *_initialize_rb_system._rb_con_ptr->get_Fq(q_f);
        temp_vec->scale( _initialize_rb_system._rb_con_ptr->get_rb_theta_expansion().eval_F_theta(q_f, mu) );
        _initialize_rb_system._rb_con_ptr->rhs->add(*temp_vec);
//        _console << *_initialize_rb_system._rb_con_ptr->rhs << std::endl;
      }
  }

  _initialize_rb_system._rb_con_ptr->matrix->close();
  _initialize_rb_system._rb_con_ptr->rhs->close();
}

void
DwarfElephantOfflineStage::setInnerProductMatrix()
{
//  if (_initialize_rb_system._qa > 1)
//  {
    for(std::set<SubdomainID>::const_iterator it = _subdomain_ids.begin();
            it != _subdomain_ids.end(); it++)
    {
      _cache_stiffness_matrix->setCachedSubdomainStiffnessMatrixContributions(*_initialize_rb_system._jacobian_subdomain[*it], *it);
      _initialize_rb_system._jacobian_subdomain[*it] ->close();
    }
//   }

//   else if (_initialize_rb_system._qa==1)
//   {
//     _initialize_rb_system._jacobian_subdomain[_initialize_rb_system._qa-1] -> close();
//     _initialize_rb_system._jacobian_subdomain[_initialize_rb_system._qa-1]->add(1, *_sys.matrix);
//   }
//
//    _cache_stiffness_matrix->setCachedStiffnessMatrixContributions(*_initialize_rb_system._inner_product_matrix);
    _initialize_rb_system._inner_product_matrix -> close();
    _initialize_rb_system._inner_product_matrix->add(1., *_sys.matrix);
}

void
DwarfElephantOfflineStage::transferAffineVectors()
{
    // Transfer the vectors
    // Transfer the data for the F vectors.
    for(unsigned int _q=0; _q<_initialize_rb_system._qf; _q++)
    {
//        _initialize_rb_system._residuals[_q]->operator=(_sys.get_vector(_residual_name));
//        _cache_stiffness_matrix->setCachedResidual(*_initialize_rb_system._residuals[_q]);
      _initialize_rb_system._residuals[_q]->close();
    }
    // Transfer the data for the output vectors.
    if (_compliant)
    {
        for(unsigned int _q=0; _q<_initialize_rb_system._ql; _q++)
        {
//          _initialize_rb_system._outputs[_q]->operator=(_sys.get_vector(_residual_name));
//          _cache_stiffness_matrix->setCachedResidual(*_initialize_rb_system._outputs[_q]);
          _initialize_rb_system._outputs[_q]->close();
        }
    }
    else if (!_compliant)
        mooseError("Currently, the code handles the compliant case, only.");
}

void
DwarfElephantOfflineStage::offlineStage()
{
    // This method performs the offline stage of the RB problem.

    // Computation of the reduced basis space.
    trainReducedBasis();
//    _initialize_rb_system._rb_con_ptr->train_reduced_basis();

    // Write the offline data to file (xdr format).
    _initialize_rb_system._rb_con_ptr->get_rb_evaluation().legacy_write_offline_data_to_files();

    // If desired, store the basis functions (xdr format).
    if (_store_basis_functions)
    {
        _initialize_rb_system._rb_con_ptr->get_rb_evaluation().write_out_basis_functions(*_initialize_rb_system._rb_con_ptr);
    }

    _initialize_rb_system._rb_con_ptr->print_basis_function_orthogonality();
}

void
DwarfElephantOfflineStage::setOnlineParameters()
{
    for (unsigned int  _q = 0; _q != _online_mu_parameters.size(); _q++)
    {
        std::string  _mu_name = "mu_" + std::to_string(_q);
        _rb_online_mu.set_value(_mu_name, _online_mu_parameters[_q]);
    }
}

void
DwarfElephantOfflineStage::initialize()
{
    if (_residual_name != "Re_non_time" && _residual_name != "Re_time")
        mooseError ("You have choosen an invalid residual_name. Valid names are 'Re_non_time' and 'Re_time'.");
}

void
DwarfElephantOfflineStage::execute()
{
    // Build the RBEvaluation object
    // Required for both the Offline and Online stage.
    DwarfElephantRBEvaluation  _rb_eval(_mesh_ptr->comm());

    // Pass a pointer of the RBEvaluation object to the
    // RBConstruction object
    _initialize_rb_system._rb_con_ptr->set_rb_evaluation(_rb_eval);

    if(_skip_vector_assembly_in_rb_system)
        transferAffineVectors();

    if(_skip_matrix_assembly_in_rb_system)
        setInnerProductMatrix();

    _console << std::endl;
    offlineStage();
    _console << std::endl;

    if(_online_stage)
    {
    #if defined(LIBMESH_HAVE_CAPNPROTO)
    RBDataDeserialization::RBEvaluationDeserialization _rb_eval_reader(_rb_eval);
    #else
    _rb_eval.legacy_read_offline_data_from_files();
    #endif

    setOnlineParameters();
    _rb_eval.set_parameters(_rb_online_mu);

    _console << "---- Online Stage ----" << std::endl;
    _rb_eval.print_parameters();

    _rb_eval.rb_solve(_online_N);

//    for (unsigned int _q = 0; _q != _initialize_rb_system._ql; _q++)
//      _console << "Output " << std::to_string(_q) << ": value = " << _rb_eval.RB_outputs[_q]
//      << ", error bound = " << _rb_eval.RB_output_error_bounds[_q] << std::endl;

    _rb_eval.read_in_basis_functions(*_initialize_rb_system._rb_con_ptr);
    _initialize_rb_system._rb_con_ptr->load_rb_solution();
    ExodusII_IO(_mesh_ptr->getMesh()).write_equation_systems("RB_sol.e", _es);
    _initialize_rb_system._rb_con_ptr->load_basis_function(0);
    ExodusII_IO(_mesh_ptr->getMesh()).write_equation_systems("bf0.e", _es);
    }
}

void
DwarfElephantOfflineStage::finalize()
{
}
