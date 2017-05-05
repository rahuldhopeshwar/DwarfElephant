/**
 * This Kernel is implements a darcy flow problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DARCY_H
#define DARCY_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Diffusion.h"

///-------------------------------------------------------------------------
// Forward Declarations
class Darcy;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<Darcy>();

///-------------------------------------------------------------------------
class Darcy : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  Darcy(const InputParameters & parameters);

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
#endif // DARCY_H
