 /**
  * The structures are defined for an elliptic PDE with the following restrictions:
  *  1. The parametere dimension p is equal to one (an extension to p = 2 will
  *    be found in RBKernelSteadyStateP2ThetatEqualMu and an extension to
  *    p = 3 will be implemented in RBKernelSteadyStateP3ThetaEqualMu).
  * 2. Theta is equal to mu (for implementing other relationships,please
  *    follow the structure of these implementation for a general usability).
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. Aq --> stiffness matrix (parameter-independent)
  * 3. Fq --> load vector (parameter-independent)
  * 4. Output
  * 5. RBThetaExpansion
  * 6. RBAssemblyExpansion
  */

#ifndef RBSTRUCTURESP1THETAEQUALMU_H
#define RBSTRUCTURESP1THETAEQUALMU_H

#include "RBKernel.h"
#include "FEProblem.h"
#include "MooseObject.h"
#include "GeneralUserObject.h"
#include "NonlinearSystem.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem_assembly.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/getpot.h"

// rbOOmit includes
#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::ElemAssembly;
using libMesh::FEInterface;
using libMesh::FEMContext;
using libMesh::Number;
using libMesh::Point;
using libMesh::RBAssemblyExpansion;
using libMesh::RBParameters;
using libMesh::RBTheta;
using libMesh::RBThetaExpansion;
using libMesh::Real;
using libMesh::RealGradient;
using libMesh::FEBase;

///------------------------------------RB-----------------------------------

///-----------------------------------THETA---------------------------------
/**
 * Please take the name convention of this package for the mu object into
 * account to ensure a gernal useability of your class.
 *
 * NOTE: IF THE PARAMETER DIMENSION EXCEEDS ONE ADDITIONAL THETA OBJECTS
 *       ARE REQUIRED.
 */

struct ThetaA0 : RBTheta
{
  virtual Number evaluate (RBParameters & _mu)
  {
    return _mu.get_value("muTest");
  }
};

///------------------------------------A------------------------------------
/**
 * Assembles the stiffness matrix by calling the corresponding RBKernel
 * function.
 *
 * NOTE: NO EDITING REQUIRED !!!
 *       If the parameter dimension execeeds one that add addtionall A's
 *       defined analoge to A0.
 */

struct A0 : ElemAssembly, RBTheta
{
  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;
    FEBase * elem_fe = libmesh_nullptr;
    JacobianBlock * _jacobian_block_ptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealGradient> > & dphi = elem_fe->get_dphi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        for (unsigned int j=0; j != n_u_dofs; j++)
          c.get_elem_jacobian()(i,j) += 1;
  }
};

///------------------------------------F------------------------------------
/**
 * Assembles the load vector by calling the corresponding RBKernel
 * function.
 *
 * NOTE: NO EDITING REQUIRED !!!
 *       If the parameter dimension execeeds one that add addtionall F's
 *       defined analoge to F0.
 */

  struct F0 : ElemAssembly//, NonlinearSystem
{
//    F0(FEProblem & problem, const std::string & name):
//      NonlinearSystem(problem, name)
//    {}

//  unsigned int _test_number;

  virtual void interior_assembly(FEMContext & c)
  {
//    FEProblem * _problem_test_ptr;
//    NonlinearSystem * _non_sys_ptr;
//    _problem_test_ptr = MooseObject::getParam<FEProblem *>("_fe_problem");
//    _non_sys_ptr = &_problem_test_ptr->getNonlinearSystem();
//    NumericVector<Number> & _residual = _non_sys_ptr->residualVector(Moose::KT_ALL);
//    RBKernel * _rb_ptr = libmesh_nullptr;
//    c.get_elem_residual()= _rb_ptr->getResidual();
    const unsigned int u_var = 0;

    FEBase * elem_fe = libmesh_nullptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    const std::vector<std::vector<Real> > & phi = elem_fe->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.get_elem_residual()(i) += JxW[qp] * (1.*phi[i][qp]);

  }

  friend class RBOutput;
  friend class RBClasses;
};

  struct F0_test : ElemAssembly//, NonlinearSystem
{
//    F0_test(FEProblem & problem, const std::string & name_in):
//      NonlinearSystem(problem, name_in)
//    {}

  virtual void interior_assembly(FEMContext & c)
  {
    RBKernel * _rb_ptr;
    std::string parameters_filename_RB = "/home/bl1/projects/moose/Conduction/comparison_different_layer_types/threeLayers_parallel_RB.i";
    GetPot infile (parameters_filename_RB);

    unsigned int _test_number = infile("Nmax",0);


//    c.get_elem_residual() = _rb_ptr-> &_local_ke;
  }

  friend class RBOutput;
  friend class RBClasses;
};

///---------------------------------OUTPUT----------------------------------
/**
 * Assembles the output vector by calling the corresponding RBKernel
 * function.
 *
 * NOTE: EDITING REQUIRED IF SPECIFIC OUTPUT DESIRED
 */

struct O0 : ElemAssembly, RBTheta
{
  virtual void interior_assembly(FEMContext & c)
  {
//    RBKernel * _rb_ptr = libmesh_nullptr;
//    c.get_elem_residual()= _rb_ptr->getResidual();
    const unsigned int u_var = 0;

    FEBase * elem_fe = libmesh_nullptr;
    c.get_element_fe(u_var, elem_fe);

    const std::vector<Real> & JxW = elem_fe->get_JxW();

    const std::vector<std::vector<Real> > & phi = elem_fe->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_var).size();

    // Now we will build the affine operator
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.get_elem_residual()(i) += JxW[qp] * (1.*phi[i][qp]);

  }

};

///----------------------------RBTHETAEXPANSION-----------------------------
/**
 * Attaches the stiffness matrix and the theta object to a structure of the
 * type RBThetatExpansion.
 *
 * NOTE: EDITING ONLY REQUIRED IF THE PARAMETER DIMENSION EXCEEDS ONE. THEN,
 *       ATTACH ADDITIONAL THETA, A, AND F OBJECTS.
 */

struct RBP1ThetaEqualMuThetaExpansion : RBThetaExpansion
{
  RBP1ThetaEqualMuThetaExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_A_theta(&_theta_a_0);
    attach_F_theta(&_rb_theta);
    attach_output_theta(&_rb_theta); // Attach output 0 theta
//    attach_output_theta(&_rb_theta); // Attach output 1 theta
//    attach_output_theta(&_rb_theta); // Attach output 2 theta
//    attach_output_theta(&_rb_theta); // Attach output 3 theta

  }
  // Member Variables
  ThetaA0 _theta_a_0;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};

///---------------------------RBASSEMBLYEXPANSION---------------------------
/**
 * Attaches the stiffness matrix, the load vector, and the output vector to
 * a structure of the type RBAssemblyExpansion.
 *
 * NOTE: EDITING ONLY REQUIRED IF: 1. THE PARAMETER DIMENSION EXCEEDS ONE.
 *                                    THEN, ATTACH ADDITIONAL A, AND F
 *                                    OBJECTS.
 *                                 2. MORE THAN ONE OUTPUT VECTOR IS DEFINED.
 *                                    THEN, ATTACH ADITIONAL OUTPUT OBJECT.
 */
struct RBP1ThetaEqualMuAssemblyExpansion : RBAssemblyExpansion
{
  RBP1ThetaEqualMuAssemblyExpansion() //:
//    L0(0.7, 0.8, 0.7, 0.8),
//    L1(0.2, 0.3, 0.7, 0.8),
//    L2(0.2, 0.3, 0.2, 0.3),
//    L3(0.7, 0.8, 0.2, 0.3)

  {
    attach_A_assembly(&A0_assembly);
    attach_F_assembly(&F0_assembly);
    attach_output_assembly(&O0_assembly);
//    attach_output_assembly(&L0);       // Attach output 0 assembly
//    attach_output_assembly(&L1);       // Attach output 1 assembly
//    attach_output_assembly(&L2);       // Attach output 2 assembly
//    attach_output_assembly(&L3);       // Attach output 3 assembly
  }

  // Member Variables
  A0 A0_assembly;
  F0 F0_assembly;
  O0 O0_assembly;
//  OutputAssembly L0;
//  OutputAssembly L1;
//  OutputAssembly L2;
//  OutputAssembly L3;

};

///---------------------------RBASSEMBLYEXPANSION---------------------------
/**
 * Attaches the stiffness matrix, the load vector, and the output vector to
 * a structure of the type RBAssemblyExpansion.
 *
 * NOTE: EDITING ONLY REQUIRED IF: 1. THE PARAMETER DIMENSION EXCEEDS ONE.
 *                                    THEN, ATTACH ADDITIONAL A, AND F
 *                                    OBJECTS.
 *                                 2. MORE THAN ONE OUTPUT VECTOR IS DEFINED.
 *                                    THEN, ATTACH ADITIONAL OUTPUT OBJECT.
 */

struct RBP1ThetaEqualMuAssemblyExpansion_test : RBAssemblyExpansion
//struct RBP1ThetaEqualMuAssemblyExpansion_test : RBAssemblyExpansion, NonlinearSystem
{
  RBP1ThetaEqualMuAssemblyExpansion_test()// :
//  RBP1ThetaEqualMuAssemblyExpansion_test(FEProblem & _problem, const std::string name_in) :
//    NonlinearSystem(_problem, name_in)
//    L0(0.7, 0.8, 0.7, 0.8),
//    L1(0.2, 0.3, 0.7, 0.8),
//    L2(0.2, 0.3, 0.2, 0.3),
//    L3(0.7, 0.8, 0.2, 0.3)

  {
    attach_A_assembly(&A0_assembly);
    attach_F_assembly(&F0_assembly);
    attach_output_assembly(&O0_assembly);
//    attach_output_assembly(&L0);       // Attach output 0 assembly
//    attach_output_assembly(&L1);       // Attach output 1 assembly
//    attach_output_assembly(&L2);       // Attach output 2 assembly
//    attach_output_assembly(&L3);       // Attach output 3 assembly
  }

  // Member Variables
  A0 A0_assembly;
  F0 F0_assembly;
  O0 O0_assembly;
//  OutputAssembly L0;
//  OutputAssembly L1;
//  OutputAssembly L2;
//  OutputAssembly L3;

};





#endif // RBSTRUCTURESP1THETAEQUALMU_H
