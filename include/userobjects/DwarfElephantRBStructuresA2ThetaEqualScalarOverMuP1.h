///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTURESA2THETAEQUALSCALAROVERMUP1_H
#define DWARFELEPHANTRBSTRUCTURESA2THETAEQUALSCALAROVERMUP1_H

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

struct DwarfElephantThetaA2ScalarOverMu : RBTheta
{
  virtual Number evaluate (const RBParameters & _mu)
  {
    Number _scalar = 1.;

    Number _theta = _scalar/_mu.get_value("mu_2");

    return _theta;
  }
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTURESA2THETAEQUALSCALAROVERMUP1_H
