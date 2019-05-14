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
#include "DwarfElephantRBKernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBDarcy;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBDarcy>();

///This Kernel is implements a darcy flow problem using the full Finite Element solution. It is included in this package for validation purposes.
class DwarfElephantRBDarcy : public DwarfElephantRBKernel
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
  Real _permeability;
  Real _norm_value_perm;
  Real _viscosity;
  Real _norm_value_visc;
  bool _gravity_term;
  Real _fluid_density;
  RealVectorValue _gravity;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBDARCY_H
