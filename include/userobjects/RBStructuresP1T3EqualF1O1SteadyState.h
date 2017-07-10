 /**
  * The structures are defined for an elliptic PDE with the following restrictions:
  *  1. The parameter dimension p is equal to one (P1).
  *  2. The number of thetas is equal to three (T3).
  *  3. Theta is equal to mu (for implementing other relationships,please
  *     follow the structure of these implementation for a general usability)
  *     (Equal).
  *  4. The problem contains one load vector (F1) and one output (O1).
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. RBThetaExpansion
  *

///-------------------------------------------------------------------------
#ifndef RBSTRUCTURESP1T3EQUALF1O1STEADYSTATE_H
#define RBSTRUCTURESP1T3EQUALF1O1STEADYSTATE_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes (RB package)
#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"

#include "RBStructuresA00ThetaIsConstantP1.h"
#include "RBStructuresA0ThetaEqualMuP1.h"
#include "RBStructuresA1ThetaEqualMuP1.h"
#include "RBStructuresA2ThetaEqualMuP1.h"


// Forward Declarations
namespace libMesh
{
 // class RBParameters;
 // class RBTheta;
  class RBThetaExpansion;
}

///----------------------------RBTHETAEXPANSION-----------------------------
/**
 * Attaches the stiffness matrix and the theta object to a structure of the
 * type RBThetatExpansion.
 *
 */

struct RBP1T3EqualF1O1SteadyStateExpansion : RBThetaExpansion
{
  RBP1T3EqualF1O1SteadyStateExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_A_theta(&_theta_a_0);
    attach_A_theta(&_theta_a_1);
    attach_A_theta(&_theta_a_2);

    attach_F_theta(&_rb_theta);
//    std::vector<RBTheta *> _thetas = {&_theta_a_0, &_theta_a_1, &_theta_a_2};
//    attach_output_theta(_thetas);
    attach_output_theta(&_rb_theta);
  }
  // Member Variables
  ThetaA0 _theta_a_0;
  ThetaA1 _theta_a_1;
  ThetaA2 _theta_a_2;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};

///-------------------------------------------------------------------------
#endif // RBSTRUCTURESP1T3EQUALF1O1STEADYSTATE_H
