///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTRUCTURESA2THETANORMMU3P1_H
#define DWARFELEPHANTRBSTRUCTURESA2THETANORMMU3P1_H

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

struct DwarfElephantThetaA2NormMu3 : RBTheta
{
  virtual Number evaluate (const RBParameters & _mu)
  {
    return _mu.get_value("mu_2")/_mu.get_value("mu_3");
  }
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTRUCTURESA2THETANORMMU3P3_H
