/**
 * SandStone inherits from the MOOSE class Material and overrides
 * computeQpProperties.
 * Within the Class SandStone the typical rock properties of a sandstone are
 * stored.
 */

//-------------------------------------------------------------------------
#ifndef DWARFELEPHANTSANDSTONE_H
#define DWARFELEPHANTSANDSTONE_H

//---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Material.h"

//-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantSandStone;

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantSandStone>();

///SandStone inherits from the MOOSE class Material and overrides computeQpProperties. Within the Class SandStone the typical rock properties of a sandstone are stored.
class DwarfElephantSandStone : public Material
{

//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantSandStone(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:
/* Methods */
  virtual void computeQpProperties() override;

/* Attributes */
  MaterialProperty<Real> & _thermal_conductivity;
  MaterialProperty<Real> & _permeability;
  MaterialProperty<Real> & _dynamic_viscosity;
  MaterialProperty<Real> & _fluid_density;
  MaterialProperty<RealVectorValue> & _gravity;
};

//-------------------------------------------------------------------------
#endif //DWARFELEPHANTSANDSTONE_H
