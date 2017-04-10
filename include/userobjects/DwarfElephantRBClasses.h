/**
 * In this class subclasses of the RBEvaluation and
 * RBConstruction class are introduced.
 *
 * DwarfElephantRBEvaluation: requires only the definition of the lower
 * coercivity constant. The value is here specified for a three layer
 * problem.
 *
 * DwarfElephantRBConstruction: In order to construct the RB System with the
 * DwarfElephantRBEvaluation subclass the method build_rb_evaluation needs to be
 * overriden.
 *
 * NOTE: ENSURE THAT THE CLASS IS INHERITING FROM THE CORRECT RBSTRUCTURES.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBCLASSES_H
#define DWARFELEPHANTRBCLASSES_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/sparse_matrix.h"

// libMesh includes (RB package)
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"

///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "RBStructuresP1Theta3ThetaEqualMu.h"

// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;

  class EquationSystems;
  class RBConstruction;
  class RBEvaluation;
}

///-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
class DwarfElephantRBConstruction : public RBConstruction
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBConstruction (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : Parent(es, name_in, number_in)
  {}

  // Destructor
  virtual ~DwarfElephantRBConstruction () { }

  // Type of the system
  typedef DwarfElephantRBConstruction _sys_type;

  // Type of the parent
  typedef RBConstruction Parent;

  // Initialize data structure
  virtual void init_data()
  {
    u_var = this->add_variable ("u", libMesh::FIRST);

    Parent::init_data();
  }

// overwrite the function to avoid negative values underneath the square root
  void enrich_RB_space()
  {
    LOG_SCOPE("enrich_RB_space()", "RBConstruction");

    NumericVector<Number> * new_bf = NumericVector<Number>::build(this->comm()).release();
    new_bf->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    *new_bf = *solution;

    for(unsigned int index=0; index<get_rb_evaluation().get_n_basis_functions(); index++)
    {
      inner_product_matrix->vector_mult(*inner_product_storage_vector, *new_bf);

      Number scalar =
        inner_product_storage_vector->dot(get_rb_evaluation().get_basis_function(index));
      new_bf->add(-scalar, get_rb_evaluation().get_basis_function(index));
    }

    // Normalize new_bf
    inner_product_matrix->vector_mult(*inner_product_storage_vector, *new_bf);
    Number new_bf_norm = std::sqrt( std::abs(inner_product_storage_vector->dot(*new_bf)) );

    if(new_bf_norm == 0.)
    {
      new_bf->zero(); // avoid potential nan's
    }
    else
    {
      new_bf->scale(1./new_bf_norm);
    }

    // load the new basis function into the basis_functions vector.
    get_rb_evaluation().basis_functions.push_back( new_bf );
}

  unsigned int u_var;

  friend class DwarfElephantOfflineStage;
};

///------------------------DWARFELEPHANTRBEVALUATION------------------------
class DwarfElephantRBEvaluation : public RBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  DwarfElephantRBEvaluation(const libMesh::Parallel::Communicator & comm):
    RBEvaluation(comm)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }

  virtual Real get_stability_lower_bound()
  {
    const RBParameters & mu = get_parameters();
    Real min_mu_1 = std::min(mu.get_value("mu_0"), mu.get_value("mu_1"));
    return std::min(min_mu_1, mu.get_value("mu_2"));
  }

  RBP1Theta3ThetaEqualMuExpansion _rb_theta_expansion;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCLASSES_H
