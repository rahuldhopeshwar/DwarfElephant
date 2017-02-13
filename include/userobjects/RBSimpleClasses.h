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

///-------------------------------------------------------------------------
#ifndef RBSIMPLECLASSES_H
#define RBSIMPLECLASSES_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/fe_base.h"

//libMesh includes (RB package)
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"

//MOOSE includes (DwarfElephant package)
#include "RBStructuresP1_3ThetaEqualMuCompliant.h"

using libMesh::DirichletBoundary;
using libMesh::EquationSystems;
using libMesh::FEMContext;
using libMesh::RBConstruction;
using libMesh::RBEvaluation;
using libMesh::Real;
using libMesh::UniquePtr;

///---------------------------RBSIMPLEEVALUATION----------------------------
class RBSimpleEvaluation : public RBEvaluation
{

//---------------------------------PUBLIC-----------------------------------
public:
  RBSimpleEvaluation(const libMesh::Parallel::Communicator & comm):
    RBEvaluation(comm)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }

  virtual Real get_stability_lower_bound(){ return 0.05; };

  RBP1_3ThetaEqualMuThetaCompliantExpansion _rb_theta_expansion;
};

///--------------------------RBSIMPLECONSTRUCTION---------------------------
class RBSimpleConstruction : public RBConstruction
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  RBSimpleConstruction (EquationSystems & es,
                              const std::string & name_in,
                              const unsigned int number_in)
    : Parent(es, name_in, number_in),
      dirichlet_bc(UniquePtr<DirichletBoundary>())
  {}

  // Destructor
  virtual ~RBSimpleConstruction () { }

  // Type of the system
  typedef RBSimpleConstruction _sys_type;

  // Type of the parent
  typedef RBConstruction Parent;

  // Initialize data structure
 virtual void init_data()
  {
    u_var = this->add_variable ("u", libMesh::FIRST);

    // Generate a DirichletBoundary object
    dirichlet_bc = build_zero_dirichlet_boundary_object();

    // Set the Dirichet boundary IDs
    // and the Dirichlet boundary variable numbers
    dirichlet_bc->b.insert(1);
    dirichlet_bc->b.insert(3);
    dirichlet_bc->variables.push_back(u_var);

    // Attach dirichlet_bc (must do this _before_ Parent::init_data)
    get_dof_map().add_dirichlet_boundary(*dirichlet_bc);

    Parent::init_data();

    // Set the rb_assembly_expansion for this Construction object.
    // The theta expansion comes from the RBEvaluation object.
    set_rb_assembly_expansion(_rb_assembly_expansion);

    // We need to define an inner product matrix for this problem
    set_inner_product_assembly(_rb_assembly_expansion.A0_assembly);
  }

  /**
   * Pre-request all relevant element data.
   */
  virtual void init_context(FEMContext & c)
  {
    // For efficiency, we should prerequest all
    // the data we will need to build the
    // linear system before doing an element loop.
    FEBase * elem_fe = libmesh_nullptr;
    c.get_element_fe(u_var, elem_fe);

    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_dphi();

    //libMesh::out << elem_fe << std::endl;
  }

  /**
   * Variable number for u.
   */
  unsigned int u_var;

  /**
   * The object that stores the "assembly" expansion of the parameter dependent PDE,
   * i.e. the objects that define how to assemble the set of parameter-independent
   * operators in the affine expansion of the PDE.
   */
  RBP1_3ThetaEqualMuAssemblyCompliantExpansion _rb_assembly_expansion;

  /**
   * The object that defines which degrees of freedom are on a Dirichlet boundary.
   */
  UniquePtr<DirichletBoundary> dirichlet_bc;
};

///-------------------------------------------------------------------------
#endif // RBSIMPLECLASSES_H
