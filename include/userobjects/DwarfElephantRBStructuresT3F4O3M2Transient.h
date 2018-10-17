 /**
  * The structures are defined for an parabolic PDE with the following restrictions:
  *  1. The number of thetas is equal to three (T3).
  *  2. The problem contains one load vector (F1) and one output (O1).
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. RBThetaExpansion
  */
///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTUREST3F4O3M2TRANSIENT_H
#define DWARFELEPHANTRBSTRUCTUREST3F4O3M2TRANSIENT_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes (RB package)
#include "libmesh/transient_rb_theta_expansion.h"
#include "libmesh/transient_rb_assembly_expansion.h"

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


// Forward Declarations
namespace libMesh
{
  class RBTransientBThetaExpansion;
}

///The structures are defined for an parabolic PDE with the following restrictions: 1. The number of thetas is equal to three (T3). 2. The problem contains one load vector (F1) and one output (O1).
/**
 * Attaches the stiffness matrix and the theta object to a structure of the
 * type RBThetatExpansion.
 *
 */

struct DwarfElephantRBT3F4O3M2TransientExpansion : TransientRBThetaExpansion
{
  DwarfElephantRBT3F4O3M2TransientExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_M_theta(&_theta_a_00);
    attach_M_theta(&_rb_theta);

    attach_A_theta(&_theta_a_0);
    attach_A_theta(&_theta_a_1);
    attach_A_theta(&_theta_a_2);

    attach_F_theta(&_theta_a_0);
    attach_F_theta(&_theta_a_1);
    attach_F_theta(&_theta_a_2);
    attach_F_theta(&_rb_theta);

    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);
    attach_output_theta(&_rb_theta);

  }
  // Member Variables
  DwarfElephantThetaA00Constant _theta_a_00;
  DwarfElephantThetaA0NormMu1 _theta_a_0;
  DwarfElephantThetaA1NormMu1 _theta_a_1;
  DwarfElephantThetaA2NormMu1 _theta_a_2;
  RBTheta _rb_theta;          // Default RBTheta object, simply returns one.
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTUREST3F4O3M2Transient_H
