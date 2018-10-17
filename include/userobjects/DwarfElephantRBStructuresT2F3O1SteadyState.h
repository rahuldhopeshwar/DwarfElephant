 /**
  * The structures are defined for an elliptic PDE with the following restrictions:
  *  1. The number of thetas is equal to two (T2).
  *  2. The problem contains three load vectors (F3) and one output (O1).
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. RBThetaExpansion
  */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTUREST2F3O1STEADYSTATE_H
#define DWARFELEPHANTRBSTRUCTUREST2F3O1STEADYSTATE_H

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


// Forward Declarations
namespace libMesh
{
 // class RBParameters;
 // class RBTheta;
  class RBThetaExpansion;
}

///The structures are defined for an elliptic PDE with the following restrictions: 1. The number of thetas is equal to two (T2). 2. The problem contains three load vectors (F3) and one output (O1).

/**
 * Attaches the stiffness matrix and the theta object to a structure of the
 * type RBThetatExpansion.
 *
 */

struct DwarfElephantRBT2F3O1SteadyStateExpansion : RBThetaExpansion
{
  DwarfElephantRBT2F3O1SteadyStateExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_A_theta(&_rb_theta);
    attach_A_theta(&_theta_a_0);

    attach_F_theta(&_rb_theta);
    attach_F_theta(&_theta_a_0);
    attach_F_theta(&_theta_a_1);

    attach_output_theta(&_rb_theta);
  }
  // Member Variables
  DwarfElephantThetaA0EqualMu0 _theta_a_0;
  DwarfElephantThetaA1EqualMu1 _theta_a_1;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTUREST2F3O1STEADYSTATE_H
