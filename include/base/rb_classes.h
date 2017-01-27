// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic
//
//     This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef RB_CLASSES_H
#define RB_CLASSES_H

#include "libmesh/rb_construction.h"
#include "libmesh/fe_base.h"
#include "libmesh/rb_evaluation.h"

// local include
//#include "assembly.h"
#include "RBStructuresP1ThetaEqualMu.h"
#include "FEProblem.h"
#include "Moose.h"
#include "SystemBase.h"
#include "NonlinearSystem.h"
#include "RBOutput.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::DirichletBoundary;
using libMesh::EquationSystems;
using libMesh::FEMContext;
using libMesh::RBConstruction;
using libMesh::RBEvaluation;
using libMesh::Real;
using libMesh::UniquePtr;

// A simple subclass of RBEvaluation, which just needs to specify
// (a lower bound for) the coercivity constant for this problem.
// For this simple convection-diffusion problem, we can set the
// coercivity constant lower bound to 0.05.

class RBOutput;

class Test : public NonlinearSystem
{
  public:

    Test(FEProblem & problem, const std::string & name):
      NonlinearSystem(problem, name)
    {
    }

    void something(){}
};

class SimpleRBEvaluation : public RBEvaluation
{
public:

  /**
   * Constructor. Just set the theta expansion.
   */
  SimpleRBEvaluation(const libMesh::Parallel::Communicator & comm) :
    RBEvaluation(comm)
  {
    set_rb_theta_expansion(_rb_theta_expansion);
  }

  /**
   * The coercivity constant is bounded below by 0.05.
   */
  virtual Real get_stability_lower_bound() { return 0.05; }

  /**
   * The object that stores the "theta" expansion of the parameter dependent PDE,
   * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
   */
  RBP1ThetaEqualMuThetaExpansion _rb_theta_expansion;

};

// A simple subclass of Construction, which just needs to override build_rb_evaluation
// in order to build a SimpleRBEvaluation object, rather than an RBEvaluation object.
class SimpleRBConstruction : public RBConstruction

{
public:

  SimpleRBConstruction (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : Parent(es, name_in, number_in),
      dirichlet_bc(UniquePtr<DirichletBoundary>())
  {}

  /**
   * Destructor.
   */
  virtual ~SimpleRBConstruction () { }

  /**
   * The type of system.
   */
  typedef SimpleRBConstruction sys_type;

  /**
   * The type of the parent.
   */
  typedef RBConstruction Parent;

  /**
   * Initialize data structures.
   */
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
  RBP1ThetaEqualMuAssemblyExpansion _rb_assembly_expansion;

  /**
   * The object that defines which degrees of freedom are on a Dirichlet boundary.
   */
  UniquePtr<DirichletBoundary> dirichlet_bc;
};

class SimpleRBConstructionTest : public RBConstruction

{
public:

  SimpleRBConstructionTest (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : Parent(es, name_in, number_in),
      dirichlet_bc(UniquePtr<DirichletBoundary>())
    {}

//  SimpleRBConstructionTest(FEProblem & _problem, const std::string & name_in):
//    NonlinearSystem(_problem, name_in)
//  {}

  /**
   * Destructor.
   */
  virtual ~SimpleRBConstructionTest () { }

  /**
   * The type of system.
   */
  typedef SimpleRBConstructionTest sys_type;

  /**
   * The type of the parent.
   */
  typedef RBConstruction Parent;

  /**
   * Initialize data structures.
   */
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
  RBP1ThetaEqualMuAssemblyExpansion _rb_assembly_expansion;

  /**
   * The object that defines which degrees of freedom are on a Dirichlet boundary.
   */
  UniquePtr<DirichletBoundary> dirichlet_bc;
};

#endif