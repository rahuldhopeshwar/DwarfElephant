/**
 * This Kernel is implements a darcy flow problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTDARCY_H
#define DWARFELEPHANTDARCY_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantDarcy;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantDarcy>();

///-------------------------------------------------------------------------
class DwarfElephantDarcy : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantDarcy(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  const MaterialProperty<Real> & _permeability;
  const MaterialProperty<Real> & _dynamic_viscosity;
  const MaterialProperty<Real> & _fluid_density;
  const MaterialProperty<RealVectorValue> & _gravity;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTDARCY_H
