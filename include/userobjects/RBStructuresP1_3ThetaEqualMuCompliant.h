 /**
  * The structures are defined for an elliptic PDE with the following restrictions:
  *  1. The parameter dimension p is equal to one. (P1)
  *  2. Theta has three different values. (_3)
  *  3. Theta is equal to mu (for implementing other relationships,please
  *     follow the structure of these implementation for a general usability).
  *     (ThetaEqualMu)
  *  4. The problem is compliant. (Compliant)
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. Aq --> stiffness matrix (parameter-independent)
  * 3. Fq --> load vector (parameter-independent)
  * 4. Output
  * 5. RBThetaExpansion
  * 6. RBAssemblyExpansion
  */

///-------------------------------------------------------------------------
#ifndef RBSTRUCTURESP1_3THETAEQUALMUCOMPLIANT_H
#define RBSTRUCTURESP1_3THETAEQUALMUCOMPLIANT_H

///---------------------------------INCLUDES--------------------------------
#include <iostream>
#include <fstream>

// libMesh includes
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
#include "libmesh/mesh_base.h"
#include "libmesh/parallel_object.h"

#include "libmesh/exodusII_io.h"
#include "libmesh/xdr_cxx.h"

// libMesh includes (RB package)
#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"

#include "RBKernel.h"
#include "InputParameters.h"
#include "FEProblemBase.h"
#include "MooseApp.h"

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
using libMesh::MeshBase;
using libMesh::ParallelObject;

class RBKernel;

///------------------------------------RB-----------------------------------

///-----------------------------------THETA---------------------------------
/**
 * Please take the name convention of this package for the mu object into
 * account to ensure a gernal useability of your class.
 */

struct ThetaA0 : RBTheta
{
  virtual Number evaluate (RBParameters & _mu)
  {
    return _mu.get_value("mu_0");
  }
};

struct ThetaA1 : RBTheta
{
  virtual Number evaluate (RBParameters & _mu)
  {
    return _mu.get_value("mu_1");
  }
};

struct ThetaA2 : RBTheta
{
  virtual Number evaluate (RBParameters & _mu)
  {
    return _mu.get_value("mu_2");
  }
};

///------------------------------------A------------------------------------
/**
 * Assembles the stiffness matrix by calling the corresponding RBKernel
 * function.
 */

struct A0 : ElemAssembly
{
  // Assemble the Laplacian operator
  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * elem_fe = libmesh_nullptr;
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
          c.get_elem_jacobian()(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
  }
};

struct A1 : ElemAssembly
{
  // Assemble the Laplacian operator
  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * elem_fe = libmesh_nullptr;
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
          c.get_elem_jacobian()(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
  }
};

struct A2 : ElemAssembly
{
  // Assemble the Laplacian operator
  virtual void interior_assembly(FEMContext & c)
  {
    const unsigned int u_var = 0;

    FEBase * elem_fe = libmesh_nullptr;
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
          c.get_elem_jacobian()(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
  }
};
///------------------------------------F------------------------------------
/**
 * Assembles the load vector by calling the corresponding RBKernel
 * function.
 *
 */

struct F0 : ElemAssembly
{
  // Source term, 1 throughout the domain
  virtual void interior_assembly(FEMContext & c)
  {
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

///---------------------------------OUTPUT----------------------------------
/**
 * Assembles the output vector by calling the corresponding RBKernel
 * function.
 */

struct O0 : ElemAssembly
{
  // Source term, 1 throughout the domain
  virtual void interior_assembly(FEMContext & c)
  {
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
 */

struct RBP1_3ThetaEqualMuThetaCompliantExpansion : RBThetaExpansion
{
  RBP1_3ThetaEqualMuThetaCompliantExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_A_theta(&_theta_a_0);
    attach_A_theta(&_theta_a_1);
    attach_A_theta(&_theta_a_2);
    attach_F_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
  }
  // Member Variables
  ThetaA0 _theta_a_0;
  ThetaA1 _theta_a_1;
  ThetaA2 _theta_a_2;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};

///---------------------------RBASSEMBLYEXPANSION---------------------------
/**
 * Attaches the stiffness matrix, the load vector, and the output vector to
 * a structure of the type RBAssemblyExpansion.
 */
struct RBP1_3ThetaEqualMuAssemblyCompliantExpansion : RBAssemblyExpansion
{
  RBP1_3ThetaEqualMuAssemblyCompliantExpansion()

  {
    attach_A_assembly(&A0_assembly);
    attach_A_assembly(&A1_assembly);
    attach_A_assembly(&A2_assembly);
    attach_F_assembly(&F0_assembly);
    attach_output_assembly(&O0_assembly);
  }

  // Member Variables
  A0 A0_assembly;
  A1 A1_assembly;
  A2 A2_assembly;
  F0 F0_assembly;
  O0 O0_assembly;
};

///-------------------------------------------------------------------------
#endif // RBSTRUCTURESP1_3THETAEQUALMUCOMPLIANT_H
