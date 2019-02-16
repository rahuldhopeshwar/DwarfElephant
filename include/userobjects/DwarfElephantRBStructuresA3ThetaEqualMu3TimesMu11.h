///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTURESA3THETAEQUALMU3TIMESMU11_H
#define DWARFELEPHANTRBSTRUCTURESA3THETAEQUALMU3TIMESMU11_H

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

/**
 * Please take the name convention of this package for the mu object into
 * account to ensure a gernal useability of your class.
 */

struct DwarfElephantThetaA3EqualMu3TimesMu11 : RBTheta
{
  virtual Number evaluate (const RBParameters & _mu)
  {
    return (_mu.get_value("mu_3")*_mu.get_value("mu_11"));
  }
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTURESA3THETAEQUALMU3TIMESMU11_H
