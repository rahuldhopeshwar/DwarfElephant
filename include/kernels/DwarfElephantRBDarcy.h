/**
 * This Kernel is implements a darcy flow problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBDARCY_H
#define DWARFELEPHANTRBDARCY_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "RBKernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBDarcy;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBDarcy>();

///-------------------------------------------------------------------------
class DwarfElephantRBDarcy : public RBKernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBDarcy(const InputParameters & parameters);

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
#endif // DWARFELEPHANTRBDARCY_H
