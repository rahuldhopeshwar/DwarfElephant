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
}

Real
DwarfElephantRBConstruction4DVar::train_reduced_basis(const bool resize_rb_eval_data)
{
  Real value = Parent::train_reduced_basis(resize_rb_eval_data);
  add_adjoint_solution(0);
  get_adjoint_solution(0).add(*solution);
  get_adjoint_solution(0).pointwise_mult(get_adjoint_solution(0),*get_output_vector(0,0));

  adjoint_solve();
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
