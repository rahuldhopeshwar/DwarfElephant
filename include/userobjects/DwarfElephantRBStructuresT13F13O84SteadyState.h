 /**
  * The structures are defined for an elliptic PDE with the following restrictions:
  *  1. The number of thetas is equal to 13 (T13).
  *  2. The problem contains 13 load vectors (F13) and 84 outputs (O84).
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. RBThetaExpansion
  */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTUREST13F13O84STEADYSTATE_H
#define DWARFELEPHANTRBSTRUCTUREST13F13O84STEADYSTATE_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes (RB package)
#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"

#include "DwarfElephantRBStructuresA00ThetaConstant.h"
#include "DwarfElephantRBStructuresA0ThetaEqualMu0.h"
#include "DwarfElephantRBStructuresA0ThetaNormMu0.h"
#include "DwarfElephantRBStructuresA0ThetaNormMu1.h"
#include "DwarfElephantRBStructuresA1ThetaEqualMu1.h"
#include "DwarfElephantRBStructuresA1ThetaNormMu0.h"
#include "DwarfElephantRBStructuresA1ThetaNormMu1.h"
#include "DwarfElephantRBStructuresA2ThetaEqualMu2.h"
#include "DwarfElephantRBStructuresA2ThetaNormMu0.h"
#include "DwarfElephantRBStructuresA2ThetaNormMu1.h"
#include "DwarfElephantRBStructuresA3ThetaEqualMu3.h"
#include "DwarfElephantRBStructuresA4ThetaEqualMu4.h"
#include "DwarfElephantRBStructuresA5ThetaEqualMu5.h"
#include "DwarfElephantRBStructuresA6ThetaEqualMu6.h"
#include "DwarfElephantRBStructuresA7ThetaEqualMu7.h"
#include "DwarfElephantRBStructuresA8ThetaEqualMu8.h"
#include "DwarfElephantRBStructuresA9ThetaEqualMu9.h"
#include "DwarfElephantRBStructuresA10ThetaEqualMu10.h"
#include "DwarfElephantRBStructuresA11ThetaEqualMu11.h"
#include "DwarfElephantRBStructuresA12ThetaEqualMu12.h"
#include "DwarfElephantRBStructuresA0ThetaEqualMu0TimesMu12.h"
#include "DwarfElephantRBStructuresA1ThetaEqualMu1TimesMu12.h"
#include "DwarfElephantRBStructuresA2ThetaEqualMu2TimesMu12.h"
#include "DwarfElephantRBStructuresA3ThetaEqualMu3TimesMu12.h"
#include "DwarfElephantRBStructuresA4ThetaEqualMu4TimesMu12.h"
#include "DwarfElephantRBStructuresA5ThetaEqualMu5TimesMu12.h"
#include "DwarfElephantRBStructuresA6ThetaEqualMu6TimesMu12.h"
#include "DwarfElephantRBStructuresA7ThetaEqualMu7TimesMu12.h"
#include "DwarfElephantRBStructuresA8ThetaEqualMu8TimesMu12.h"
#include "DwarfElephantRBStructuresA9ThetaEqualMu9TimesMu12.h"
#include "DwarfElephantRBStructuresA10ThetaEqualMu10TimesMu12.h"
#include "DwarfElephantRBStructuresA11ThetaEqualMu11TimesMu12.h"


// Forward Declarations
namespace libMesh
{
 // class RBParameters;
 // class RBTheta;
  class RBThetaExpansion;
}

///The structures are defined for an elliptic PDE with the following restrictions: 1. The number of thetas is equal to nine (T9). 2. The problem contains two load vectors (F2) and 80 outputs (O80).

/**
 * Attaches the stiffness matrix and the theta object to a structure of the
 * type RBThetatExpansion.
 *
 */

struct DwarfElephantRBT13F13O84SteadyStateExpansion : RBThetaExpansion
{
  DwarfElephantRBT13F13O84SteadyStateExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_A_theta(&_theta_a_0);
    attach_A_theta(&_theta_a_1);
    attach_A_theta(&_theta_a_2);
    attach_A_theta(&_theta_a_3);
    attach_A_theta(&_theta_a_4);
    attach_A_theta(&_theta_a_5);
    attach_A_theta(&_theta_a_6);
    attach_A_theta(&_theta_a_7);
    attach_A_theta(&_theta_a_8);
    attach_A_theta(&_theta_a_9);
    attach_A_theta(&_theta_a_10);
    attach_A_theta(&_theta_a_11);

    attach_F_theta(&_theta_a_12);
    attach_F_theta(&_theta_f_0);
    attach_F_theta(&_theta_f_1);
    attach_F_theta(&_theta_f_2);
    attach_F_theta(&_theta_f_3);
    attach_F_theta(&_theta_f_4);
    attach_F_theta(&_theta_f_5);
    attach_F_theta(&_theta_f_6);
    attach_F_theta(&_theta_f_7);
    attach_F_theta(&_theta_f_8);
    attach_F_theta(&_theta_f_9);
    attach_F_theta(&_theta_f_10);
    attach_F_theta(&_theta_f_11);

    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
  }
  // Member Variables
  DwarfElephantThetaA0EqualMu0 _theta_a_0;
  DwarfElephantThetaA1EqualMu1 _theta_a_1;
  DwarfElephantThetaA2EqualMu2 _theta_a_2;
  DwarfElephantThetaA3EqualMu3 _theta_a_3;
  DwarfElephantThetaA4EqualMu4 _theta_a_4;
  DwarfElephantThetaA5EqualMu5 _theta_a_5;
  DwarfElephantThetaA6EqualMu6 _theta_a_6;
  DwarfElephantThetaA7EqualMu7 _theta_a_7;
  DwarfElephantThetaA8EqualMu8 _theta_a_8;
  DwarfElephantThetaA9EqualMu9 _theta_a_9;
  DwarfElephantThetaA10EqualMu10 _theta_a_10;
  DwarfElephantThetaA11EqualMu11 _theta_a_11;
  DwarfElephantThetaA12EqualMu12 _theta_a_12;
  DwarfElephantThetaA0EqualMu0TimesMu12 _theta_f_0;
  DwarfElephantThetaA1EqualMu1TimesMu12 _theta_f_1;
  DwarfElephantThetaA2EqualMu2TimesMu12 _theta_f_2;
  DwarfElephantThetaA3EqualMu3TimesMu12 _theta_f_3;
  DwarfElephantThetaA4EqualMu4TimesMu12 _theta_f_4;
  DwarfElephantThetaA5EqualMu5TimesMu12 _theta_f_5;
  DwarfElephantThetaA6EqualMu6TimesMu12 _theta_f_6;
  DwarfElephantThetaA7EqualMu7TimesMu12 _theta_f_7;
  DwarfElephantThetaA8EqualMu8TimesMu12 _theta_f_8;
  DwarfElephantThetaA9EqualMu9TimesMu12 _theta_f_9;
  DwarfElephantThetaA10EqualMu10TimesMu12 _theta_f_10;
  DwarfElephantThetaA11EqualMu11TimesMu12 _theta_f_11;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTUREST13F13O84STEADYSTATE_H
