 /**
  * The structures are defined for an elliptic PDE with the following restrictions:
  *  1. The parameter dimension p is equal to one. (P1)
  *  2. Theta has three different values. (_3)
  *  3. Theta is equal to mu (for implementing other relationships,please
  *     follow the structure of these implementation for a general usability).
  *     (ThetaEqualMu)
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. RBThetaExpansion
  */

///-------------------------------------------------------------------------
#ifndef RBSTRUCTURESP1THETA5THETAEQUALMUTRANSIENT_H
#define RBSTRUCTURESP1THETA5THETAEQUALMUTRANSIENT_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes (RB package)
#include "libmesh/rb_theta.h"
#include "libmesh/transient_rb_theta_expansion.h"
//#include "libmesh/rb_assembly_expansion.h"

#include "RBStructuresA0ThetaEqualMuP1.h"
#include "RBStructuresA1ThetaEqualMuP1.h"
#include "RBStructuresA2ThetaEqualMuP1.h"
#include "RBStructuresA3ThetaEqualMuP1.h"
#include "RBStructuresA4ThetaEqualMuP1.h"


// Forward Declarations
namespace libMesh
{
 // class RBParameters;
 // class RBTheta;
  class TransientRBThetaExpansion;
}

///----------------------------RBTHETAEXPANSION-----------------------------
/**
 * Attaches the stiffness matrix and the theta object to a structure of the
 * type RBThetatExpansion.
 *
 */

struct RBP1Theta5ThetaEqualMuExpansionTransient : TransientRBThetaExpansion
{
  RBP1Theta5ThetaEqualMuExpansionTransient()
  {
    // Setting up the RBThetaExpansion object
    attach_A_theta(&_theta_a_0);
    attach_A_theta(&_theta_a_1);
    attach_A_theta(&_theta_a_2);
    attach_A_theta(&_theta_a_3);
    attach_A_theta(&_theta_a_4);

    attach_F_theta(&_theta_a_0);
    attach_F_theta(&_theta_a_1);
    attach_F_theta(&_theta_a_2);
    attach_F_theta(&_theta_a_3);
    attach_F_theta(&_theta_a_4);
//    attach_F_theta(&_rb_theta);
//    attach_output_theta(&_rb_theta);

    std::vector <RBTheta *> _thetas = {&_rb_theta, &_rb_theta, &_rb_theta, &_rb_theta, &_rb_theta};
    attach_output_theta(_thetas);
//    attach_output_theta(&_theta_a_0);
//    attach_output_theta(&_theta_a_1);
//    attach_output_theta(&_theta_a_2);
  }
  // Member Variables
  ThetaA0 _theta_a_0;
  ThetaA1 _theta_a_1;
  ThetaA2 _theta_a_2;
  ThetaA3 _theta_a_3;
  ThetaA4 _theta_a_4;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};

///-------------------------------------------------------------------------
#endif // RBSTRUCTURESP1THETA5THETAEQUALMUTRANSIENT_H
