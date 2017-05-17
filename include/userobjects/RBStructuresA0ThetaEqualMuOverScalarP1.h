///-------------------------------------------------------------------------
#ifndef RBSTRUCTURESA0THETAEQUALMUOVERSCALARP1_H
#define RBSTRUCTURESA0THETAEQUALMUOVERSCALARP1_H

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

struct ThetaA0 : RBTheta
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
#endif // RBSTRUCTURESA0THETAEQUALMUOVERSCALARP1_H
