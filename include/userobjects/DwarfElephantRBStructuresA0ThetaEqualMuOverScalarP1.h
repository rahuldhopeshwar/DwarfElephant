///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTURESA0THETAEQUALMUOVERSCALARP1_H
#define DWARFELEPHANTRBSTRUCTURESA0THETAEQUALMUOVERSCALARP1_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes (RB package)
#include "libmesh/rb_theta.h"
#include "libmesh/rb_assembly_expansion.h"


// Forward Declarations
namespace libMesh
{
  class RBParameters;
  class RBTheta;
}

///-----------------------------------THETA---------------------------------
/**
 * Please take the name convention of this package for the mu object into
 * account to ensure a gernal useability of your class.
 */

struct DwarfElephantThetaA0OverScalar : RBTheta
{
  virtual Number evaluate (const RBParameters & _mu)
  {
    Number _mu = _mu.get_value("mu_0");
    Number _scalar = 0.001145;

    Number _theta = _mu/_scalar;

    return _theta;
  }
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTURESA0THETAEQUALMUOVERSCALARP1_H
