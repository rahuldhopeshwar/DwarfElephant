/**
 * Fault inherits from the MOOSE class Material and overrides
 * computeQpProperties.
 * Within the Class Fault the typical properties of a fault are
 * stored.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFAULT_H
#define DWARFELEPHANTFAULT_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Material.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFault;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFault>();

///-------------------------------------------------------------------------
class DwarfElephantFault : public Material
{

//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFault(const InputParameters & parameters);

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

///-------------------------------------------------------------------------
#endif //DWARFELEPHANTFAULT_H
