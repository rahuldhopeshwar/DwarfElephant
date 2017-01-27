#ifndef SHALE_H_
#define SHALE_H_

//---------------------------------INCLUDE---------------------------------
// MOOSE includes
#include "Material.h"

//-------------------------------------------------------------------------
class Shale;

//-------------------------------------------------------------------------
template<>
InputParameters validParams<Shale>();

//-------------------------------------------------------------------------
/**
 * Shale inherits from the MOOSE class  Material and overrides computeQpProperties.
 * Within the class Shale the typical rock properties of a shale are stored.
 */

//-------------------------------------------------------------------------
class Shale : public Material
{

//---------------------------------PUBLIC----------------------------------
public:
  Shale(const InputParameters & parameters);

//-------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual void computeQpProperties() override;

  /* Attributes */
  /// The thermal conductivity (lambda)
  MaterialProperty<Real> & _lambda;

  // RB parameters for the thermal conductivity
  MaterialProperty<RealVectorValue> & _rb_lambda;
  MaterialProperty<Real> & _online_lambda;
};

//-------------------------------------------------------------------------
#endif //SHALE_H
