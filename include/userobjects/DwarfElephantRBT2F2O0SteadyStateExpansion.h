# pragma once

#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"

#include "DwarfElephantRBStructuresA0ThetaEqualMu0.h"
#include "DwarfElephantRBStructuresA1ThetaEqualMu1.h"
#include "DwarfElephantRBStructuresA00ThetaConstant.h"

namespace libMesh
{
  class RBThetaExpansion;
}

struct DwarfElephantThetaA00UnitConstant : DwarfElephantThetaA00Constant
{
  Number evaluate(const RBParameters &)
  {
    return 1.0;
  }
};

struct DwarfElephantRBT2F2O0SteadyStateExpansion : RBThetaExpansion
{
  DwarfElephantRBT2F2O0SteadyStateExpansion()
  {
    // Setting up the RBThetaExpansion object
    attach_A_theta(&_theta_a_0);
    attach_A_theta(&_theta_a_1);
    // attach_A_theta(&_theta_a_lm1);
    // attach_A_theta(&_theta_a_lm2);

    attach_F_theta(&_rb_theta);
    // attach_F_theta(&_rb_theta);
  }
  // Member Variables
  DwarfElephantThetaA0EqualMu0 _theta_a_0;
  DwarfElephantThetaA1EqualMu1 _theta_a_1;
  DwarfElephantThetaA00UnitConstant _theta_a_lm1;
  DwarfElephantThetaA00UnitConstant _theta_a_lm2;
  RBTheta _rb_theta;         // Default RBTheta object, simply returns one.
};
