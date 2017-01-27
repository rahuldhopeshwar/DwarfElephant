#ifndef SANDSTONE_H_
#define SANDSTONE_H_

//---------------------------------INCLUDE---------------------------------
// MOOSE includes
#include "Material.h"

//-------------------------------------------------------------------------
class SandStone;

//-------------------------------------------------------------------------
template<>

//-------------------------------------------------------------------------
InputParameters validParams<SandStone>();

//-------------------------------------------------------------------------
/**
 * SandStone inherits from the MOOSE class Material and overrides computeQpProperties.
 * Within the Class SandStone the typical rock properties of a sandstone are stored.
 */

//-------------------------------------------------------------------------
class SandStone : public Material
{

//---------------------------------PUBLIC----------------------------------
public:
  SandStone(const InputParameters & parameters);

//-------------------------------PROTECTED---------------------------------
//protected:

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
#endif //SANDSTONE_H
