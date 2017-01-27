#ifndef RBKERNEL_H
#define RBKERNEL_H

/**
 * This Kernel is required for the RB approach. It does not include any
 * RB calculations but it extracts the necessary Jacobian matrix and the
 * residuals.

 * Note that at the moment the RB approach in MOOSE is developed for an
 * elliptic PDE.
 **/

///---------------------------------INCLUDE---------------------------------

// MOOSE includes
#include "Kernel.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem_assembly.h"
#include "libmesh/quadrature_gauss.h"

// rbOOmit includes
#include "libmesh/rb_theta.h"
#include "libmesh/rb_theta_expansion.h"
#include "libmesh/rb_assembly_expansion.h"
#include "libmesh/rb_parameters.h"

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

///---------------------------FOWARD DECLARTIONS----------------------------

class RBKernel;
///-------------------------------------------------------------------------
template<>

///----------------------------INPUT PARAMETERS-----------------------------
InputParameters validParams<RBKernel>();

///------------------------------KERNEL CLASS-------------------------------
class RBKernel : public Kernel
{
///---------------------------------PUBLIC----------------------------------
public:
  RBKernel(const InputParameters & parameters);

///-------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  // Calculation of the required PDEs
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

//  const MaterialProperty<RealVectorValue> & getMu0();

  const MaterialProperty<RealVectorValue> &_mu0;


public:
  unsigned int _muTest;
  unsigned int _test;
  DenseMatrix<Number> & getJacobian();
  DenseVector<Number> & getResidual();

//  Real getMu0();

///------------------------------------RB-----------------------------------
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

struct A0 : ElemAssembly
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

    RBKernel * _rb_ptr;
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
struct F0 : ElemAssembly
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

///---------------------------------OUTPUT----------------------------------
/**
 * Assembles the output vector by calling the corresponding RBKernel
 * function.
 *
 * NOTE: EDITING REQUIRED IF SPECIFIC OUTPUT DESIRED
 */

//struct OutputAssembly : ElemAssembly
struct O0 : ElemAssembly
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


///----------------------------------------------------------------------------
friend class RBOutput;
friend class RBConductionBase;
friend struct test;
//friend struct A0;
//friend struct F0;
//friend struct O0;
};
#endif // RBKERNEL_H
