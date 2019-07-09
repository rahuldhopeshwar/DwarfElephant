 /**
  * The structures are defined for an elliptic PDE with the following restrictions:
  *  1. The number of thetas is equal to two (T2).
  *  2. The problem contains one load vector (F1) and five output (O5).
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. RBThetaExpansion
  */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTUREST2F1O5STEADYSTATE_H
#define DWARFELEPHANTRBSTRUCTUREST2F1O5STEADYSTATE_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes (RB package)
#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"

#include "DwarfElephantRBStructuresA0ThetaEqualScalarDividedByMu0.h"
#include "DwarfElephantRBStructuresA1ThetaEqualScalarDividedByMu1.h"


// Forward Declarations
namespace libMesh
{
 // class RBParameters;
 // class RBTheta;
  class RBThetaExpansion;
}

///The structures are defined for an elliptic PDE with the following restrictions: 1. The number of thetas is equal to two (T2). 2. The problem contains one load vector (F1) and one output (O1).

/**
 * Attaches the stiffness matrix and the theta object to a structure of the
 * type RBThetatExpansion.
 *
 */

struct DwarfElephantRBT2F1O5SteadyStateExpansion : RBThetaExpansion
{
  DwarfElephantRBT2F1O5SteadyStateExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_A_theta(&_theta_a_0);
    attach_A_theta(&_theta_a_1);

    attach_F_theta(&_rb_theta);

    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);

  }
  // Member Variables
  DwarfElephantThetaA0EqualScalarDividedByMu0 _theta_a_0;
  DwarfElephantThetaA1EqualScalarDividedByMu1 _theta_a_1;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTUREST2F1O5STEADYSTATE_H
