#ifndef DWARFELEPHANTEIMSTRUCTURES_H
#define DWARFELEPHANTEIMSTRUCTURES_H

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/dof_map.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"

// rbOOmit includes
#include "libmesh/rb_assembly_expansion.h"
#include "libmesh/rb_eim_theta.h"
#include "libmesh/rb_parametrized_function.h"


// MOOSE includes
#include "Function.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::ElemAssembly;
using libMesh::FEMContext;
using libMesh::Number;
using libMesh::RBAssemblyExpansion;
using libMesh::RBEIMAssembly;
using libMesh::RBEIMConstruction;
using libMesh::RBParametrizedFunction;
using libMesh::RBParameters;
using libMesh::RBTheta;
using libMesh::RBThetaExpansion;
using libMesh::Real;
using libMesh::RealGradient;
using libMesh::Elem;
using libMesh::FEBase;


struct ShiftedGaussian : public RBParametrizedFunction
{
  virtual Number evaluate(const RBParameters & mu,
                          const Point & p,
                          const Elem &)
  {
    Real center_x = mu.get_value("mu_0");
    Real center_y = mu.get_value("mu_1");
    return exp(-2.*(pow(center_x-p(0),2.) + pow(center_y-p(1),2.)));
  }
};

// Expansion of the PDE operator
struct DwarfElephantEIMThetaA0 : RBTheta { virtual Number evaluate(const RBParameters &) { return 0.05;  } };

// Define an RBThetaExpansion class for this PDE
struct DwarfElephantEIMTestRBThetaExpansion : RBThetaExpansion
{
  /**
   * Constructor.
   */
  DwarfElephantEIMTestRBThetaExpansion()
  {
    attach_A_theta(&theta_a_0);
  }

  // The RBTheta member variables
  DwarfElephantEIMThetaA0 theta_a_0;
};

struct DwarfElephantEIMFAssembly : RBEIMAssembly
{
  DwarfElephantEIMFAssembly(RBEIMConstruction & rb_eim_con_in, unsigned int basis_function_index_in) : 
    RBEIMAssembly(rb_eim_con_in, basis_function_index_in)
  {
  }

  void get_eim_basis_function_values(const Elem * _elem, const QBase * _qrule,std::vector<Number>  & eim_values)
  {
    evaluate_basis_function(0,*_elem, *_qrule,eim_values);
  }  

  virtual void interior_assembly(FEMContext & c)
  {
  }
};

#endif //DWARFELEPHANTEIMSTRUCTURES_H

