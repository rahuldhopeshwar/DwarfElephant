/**
 * Theta is a RB material class. Within this class the theta object of
 * the differnt subdomains should be defined.
 */

///-------------------------------------------------------------------------
#ifndef THETAOBJECT_H
#define THETAOBJECT_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Material.h"
#include "libmesh/rb_parameters.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class RBParameters;
}

class ThetaObject;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<ThetaObject>();

///-------------------------------------------------------------------------
class ThetaObject :
  public Material,
  public RBParameters
{

//----------------------------------PUBLIC----------------------------------
public:
  ThetaObject(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:
/* Methods */
  virtual void computeQpProperties() override;

/* Attributes */
//  RBParameters _mu;
};

///-------------------------------------------------------------------------
#endif //THETAOBJECT_H
