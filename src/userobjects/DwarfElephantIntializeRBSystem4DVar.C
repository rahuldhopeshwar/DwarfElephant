 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystem4DVar.h"

registerMooseObject("DwarfElephantApp", DwarfElephantInitializeRBSystem4DVar);

template<>
InputParameters validParams<DwarfElephantInitializeRBSystem4DVar>()
{
  InputParameters params = validParams<DwarfElephantInitializeRBSystemTransient>();
  params.addParam<unsigned int>("n_qois", "Defines the number of quantities of interest.");
  params.addParam<std::vector<Real>>("qoi_weights","A vector assiging the weigths for the qois.");

  return params;
}

DwarfElephantInitializeRBSystem4DVar::DwarfElephantInitializeRBSystem4DVar(const InputParameters & params):
  DwarfElephantInitializeRBSystemTransient(params),
  _n_qois(getParam<unsigned int>("n_qois")),
  _qoi_weights(getParam<std::vector<Real>>("qoi_weights"))
{
  if(_qoi_weights.size()!= _n_qois)
    mooseError("The number of qoi weights and qois is not consistent.");
}

void
DwarfElephantInitializeRBSystem4DVar::processParameters()
{
  DwarfElephantInitializeRBSystemTransient::processParameters();
}

// void
// DwarfElephantInitializeRBSystem4DVar::initializeOfflineStage()
// {
//   // Get and process the necessary input parameters for the
//   // offline stage
//   //  _rb_con_ptr->process_parameters_file(_parameters_filename);
//   processParameters();
//
//   // Print the system informations for the RBConstruction system.
//   _rb_con_ptr->print_info();
//
//   // Initialize the RB construction. Note, we skip the matrix and vector
//   // assembly, since this is already done by MOOSE.
//   _rb_con_ptr->initialize_rb_construction(_skip_matrix_assembly_in_rb_system, _skip_vector_assembly_in_rb_system);
//
//   // Save the A's, F's and output vectors from the RBConstruction class in pointers.
//   // This additional saving of the pointers is required because otherwise a the RBEvaluation object has
//   // to be set again in the RBKernel.
//
//   // Define size of all new parameters.
//   _jacobian_subdomain.resize(_qa);
//   _mass_matrix_subdomain.resize(_qm);
//   _residuals.resize(_qf);
//   _outputs.resize(_n_outputs);
//
//   if(_parameter_dependent_IC)
//     _inital_conditions.resize(_q_ic);
//
//   for (unsigned int i=0; i < _n_outputs; i++)
//     _outputs[i].resize(_ql[i]);
//
//   // Get the correct matrices from the RB System.
//
//   // Eliminates error message for the initialization of new non-zero entries
//   // For the future: change SparseMatrix pattern (increases efficency)
//   _inner_product_matrix = _rb_con_ptr->get_inner_product_matrix();
//   PetscMatrix<Number> * _petsc_inner_matrix = dynamic_cast<PetscMatrix<Number>* > (_inner_product_matrix);
//   MatSetOption(_petsc_inner_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//
//   // Eliminates error message for the initialization of new non-zero entries
//   // For the future: change SparseMatrix pattern (increases efficency)
//   _L2_matrix = _rb_con_ptr->L2_matrix.get();
//   PetscMatrix<Number> * _petsc_L2_matrix = dynamic_cast<PetscMatrix<Number>* > (_L2_matrix);
//   MatSetOption(_petsc_L2_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//
//   for (unsigned int _q=0; _q < _qa; _q++)
//   {
//     _jacobian_subdomain[_q] = _rb_con_ptr->get_Aq(_q);
//
//     // Eliminates error message for the initialization of new non-zero entries
//     // For the future: change SparseMatrix pattern (increases efficency)
//     PetscMatrix<Number> * _petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_jacobian_subdomain[_q]);
//     MatSetOption(_petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//   }
//
//   for (unsigned int _q=0; _q < _qm; _q++)
//   {
//     // in case you are using a libMesh version older than Dec 6, 2017 use the following
//     // line
//     // _mass_matrix_subdomain[_q] = _rb_con_ptr->M_q_vector[_q];
//     _mass_matrix_subdomain[_q] = _rb_con_ptr->get_M_q(_q);
//
//     // Eliminates error message for the initialization of new non-zero entries
//     // For the future: change SparseMatrix pattern (increases efficency)
//     PetscMatrix<Number> * _petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (_mass_matrix_subdomain[_q]);
//     MatSetOption(_petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//    }
//
//    // Get the correct vectors from the RB System.
//    for (unsigned int _q=0; _q < _qf; _q++)
//      _residuals[_q] = _rb_con_ptr->get_Fq(_q);
//
//   if(_parameter_dependent_IC)
//   {
//     DwarfElephantRBConstructionTransient * _dwarf_elephant_rb_con_ptr = dynamic_cast<DwarfElephantRBConstructionTransient * > (_rb_con_ptr);
//     for (unsigned int _q=0; _q < _q_ic; _q++)
//       _inital_conditions[_q] = _dwarf_elephant_rb_con_ptr->get_IC_q(_q);
//   }
//
//    for (unsigned int i=0; i < _n_outputs; i++)
//      for (unsigned int _q=0; _q < _ql[i]; _q++)
//        _outputs[i][_q] = _rb_con_ptr->get_output_vector(i,_q);
// }

// void
// DwarfElephantInitializeRBSystem4DVar::initialize()
// {
// }

void
DwarfElephantInitializeRBSystem4DVar::execute()
{
  // Define the parameter file for the libMesh functions.
  // GetPot infile (_parameters_filename);

  // Add a new equation system for the RB construction.
  _rb_con_ptr = &_es.add_system<DwarfElephantRBConstruction4DVar> ("RBSystem");

  if (_offline_stage){
    // Intialization of the added equation system
    _rb_con_ptr->init();
    _es.update();

    DwarfElephantRBEvaluation4DVar _rb_eval(_mesh_ptr->comm(), _fe_problem);
    // Pass a pointer of the RBEvaluation object to the
    // RBConstruction object
    _rb_con_ptr->set_rb_evaluation(_rb_eval);

    TransientRBThetaExpansion & _trans_theta_expansion = cast_ref<TransientRBThetaExpansion &>(_rb_con_ptr->get_rb_theta_expansion());

    // Get number of attached parameters.
    _n_outputs = _rb_con_ptr->get_rb_theta_expansion().get_n_outputs();
    _ql.resize(_n_outputs);
    _qa = _rb_con_ptr->get_rb_theta_expansion().get_n_A_terms();
    _qm = _trans_theta_expansion.get_n_M_terms();
    _qf = _rb_con_ptr->get_rb_theta_expansion().get_n_F_terms();

    if(_parameter_dependent_IC)
    {
      DwarfElephantRBTransientThetaExpansion & dwarf_elephant_trans_theta_expansion =
        cast_ref<DwarfElephantRBTransientThetaExpansion &>(_rb_con_ptr->get_rb_theta_expansion());

        _q_ic = dwarf_elephant_trans_theta_expansion.get_n_IC_terms();
    }

    for(unsigned int i=0; i < _n_outputs; i++)
    _ql[i] = _rb_con_ptr->get_rb_theta_expansion().get_n_output_terms(i);

    // Initialize required matrices and vectors.
      initializeOfflineStage();
  }
}

void
DwarfElephantInitializeRBSystem4DVar::finalize()
{
  // // Way to define QoI for the adjoint system
  // // Declare a QoISet object, we need this object to set weights for our QoI error contributions
  // QoISet qois;
  // _rb_con_ptr->qoi.resize(_n_qois);
  //
  // std::vector<unsigned int> qoi_indices;
  // for (unsigned int i=0; i<_n_qois; i++)
  //   qoi_indices.push_back(i);
  //
  // qois.add_indices(qoi_indices);
  //
  // for (unsigned int i=0; i<_n_qois; i++)
  //   qois.set_weight(i, _qoi_weights[i]);
}

// std::vector<std::vector<NumericVector <Number> *> >
// DwarfElephantInitializeRBSystemTransient::getOutputs() const
// {
//   return _outputs;
// }
