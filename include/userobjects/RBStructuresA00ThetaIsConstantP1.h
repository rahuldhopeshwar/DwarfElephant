///-------------------------------------------------------------------------
#ifndef RBSTRUCTURESA00THETAISCONSTANTP1_H
#define RBSTRUCTURESA00THETAISCONSTANTP1_H

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

struct ThetaA00 : RBTheta
{
  virtual Number evaluate (const RBParameters & )
  {
    Real _constant = 0.05;

    return _constant;
  }
};

///-------------------------------------------------------------------------
#endif // RBSTRUCTURESA00THETAISCONSTANTP1_H
