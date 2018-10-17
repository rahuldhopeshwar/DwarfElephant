 /**
  * The structures are defined for an parabolic PDE with the following restrictions:.
  *  1. The number of thetas is equal to one (T1).
  *  2. The problem contains one load vector (F1) and one output (O1).
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. RBThetaExpansion
  */
///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTUREST1F1O1M1TRANSIENT_H
#define DWARFELEPHANTRBSTRUCTUREST1F1O1M1TRANSIENT_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes (RB package)
#include "libmesh/transient_rb_theta_expansion.h"
#include "libmesh/transient_rb_assembly_expansion.h"

#include "DwarfElephantRBStructuresA00ThetaConstant.h"
#include "DwarfElephantRBStructuresA0ThetaEqualMu0.h"


// Forward Declarations
namespace libMesh
{
 // class RBParameters;
 // class RBTheta;
  class RBTransientBThetaExpansion;
}

/**
 * Attaches the stiffness matrix and the theta object to a structure of the
 * type RBThetatExpansion.
 *
 */

 ///The structures are defined for an parabolic PDE with the following restrictions: 1. The number of thetas is equal to one (T1). 2. The problem contains one load vector (F1) and one output (O1).
struct DwarfElephantRBT1F1O1M1TransientExpansion : TransientRBThetaExpansion
{
  DwarfElephantRBT1F1O1M1TransientExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_M_theta(&_rb_theta);

    attach_A_theta(&_theta_a_0);

    attach_F_theta(&_rb_theta);

    attach_output_theta(&_rb_theta);

  }
  // Member Variables
  DwarfElephantThetaA0EqualMu0 _theta_a_0;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTUREST1F1O1M1Transient_H
