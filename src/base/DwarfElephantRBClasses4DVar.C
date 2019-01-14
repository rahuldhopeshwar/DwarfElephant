/**
 * In this class simplified subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
 *
 * NOTE: ENSURE THAT THE CLASS IS INHERITING FROM THE CORRECT RBSTRUCTURES.
 */

//-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
#include "DwarfElephantRBClasses4DVar.h"

DwarfElephantRBConstruction4DVar::DwarfElephantRBConstruction4DVar (EquationSystems & es,
                      const std::string & name_in,
                      const unsigned int number_in)
  : Parent(es, name_in, number_in)
  // n_qois(0),
  // qoi_weights(0)
{}

void
DwarfElephantRBConstruction4DVar::init_data()
{
  Parent::init_data();
  // Cq_vector.clear();
}

void DwarfElephantRBConstruction4DVar::allocate_data_structures()
{
  Parent::allocate_data_structures();

  // Cq_vector.resize(get_rb_theta_expansion().get_n_outputs());
  //
  // for (unsigned int q=0; q<get_rb_theta_expansion().get_n_outputs(); q++)
  //   {
  //     // Initialize the memory for the matrices
  //     Cq_vector[q] = NumericVector<Number>::build(this->comm());
  //     Cq_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
  //   }
}

Real
DwarfElephantRBConstruction4DVar::train_reduced_basis(const bool resize_rb_eval_data)
{
  Real value = Parent::train_reduced_basis(resize_rb_eval_data);
  add_adjoint_solution(0);
  // get_adjoint_solution(0).add(*solution);
  // get_adjoint_solution(0).pointwise_mult(get_adjoint_solution(0),*get_output_vector(0,0));
  adjoint_solve();

  // for (unsigned int q=0; q<get_rb_theta_expansion().get_n_outputs(); q++)
  // {
  //   Cq_vector[q]->add(*get_output_vector(q,0));
  // }
  //
  // Cq_vector*get_output_vector(0,0);

  // Construct the rhs for the adjoints
  // Note: Currently, we are assuming that the covariance matrix D is equal to
  // the identity matrix. Meaning that all measurement errors are uncorrelated.
  // TODO: extend approach for correlated measurement errors --> introduce D
  DenseVector<Number> temp_vec;
  temp_vec.resize(get_rb_theta_expansion().get_n_outputs());
  for (unsigned int q=0; q<get_rb_theta_expansion().get_n_outputs(); q++)
  {
    temp_vec(q) = get_output_vector(q,0)->dot(*get_output_vector(q,0)); //replace second output vector with state
  }

  temp_vec.scale(-1);
  temp_vec-=obs_data_all_k[0];

  return value;
}

// void
// DwarfElephantRBConstruction4DVar::set_n_qois(unsigned int n_qois_in)
// {
//   this->n_qois = n_qois_in;
// }
//
// void
// DwarfElephantRBConstruction4DVar::set_qoi_weights(std::vector<Real> qoi_weights_in)
// {
//   this->qoi_weights = qoi_weights_in;
// }

//------------------------DWARFELEPHANTRBEVALUATION------------------------
#include "libmesh/xdr_cxx.h"

  DwarfElephantRBEvaluation4DVar::DwarfElephantRBEvaluation4DVar(const libMesh::Parallel::Communicator & comm, FEProblemBase & fe_problem):
    DwarfElephantRBEvaluationTransient(comm, fe_problem),
    fe_problem(fe_problem)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }
