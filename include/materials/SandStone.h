/**
 * SandStone inherits from the MOOSE class Material and overrides
 * computeQpProperties.
 * Within the Class SandStone the typical rock properties of a sandstone are
 * stored.
 */

///-------------------------------------------------------------------------
#ifndef SANDSTONE_H
#define SANDSTONE_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Material.h"

///-------------------------------------------------------------------------
// Forward Declarations
class SandStone;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<SandStone>();

///-------------------------------------------------------------------------
class SandStone : public Material
{

//----------------------------------PUBLIC----------------------------------
public:
  SandStone(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:
/* Methods */
  virtual void computeQpProperties() override;

/* Attributes */
  /// The thermal conductivity (lambda)
  MaterialProperty<Real> & _lambda;

  // RB parameters for the thermal conductivity
  MaterialProperty<RealVectorValue> & _rb_lambda;
};

///-------------------------------------------------------------------------
#endif //SANDSTONE_H
