///-------------------------------------------------------------------------
#ifndef RBSTRUCTURESA3THETAEQUALMUP1_H
#define RBSTRUCTURESA3THETAEQUALMUP1_H

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

struct ThetaA3 : RBTheta
{
  virtual Number evaluate (const RBParameters & _mu)
  {
    return _mu.get_value("mu_3");
  }
};

///-------------------------------------------------------------------------
#endif // RBSTRUCTURESA3THETAEQUALMUP1_H