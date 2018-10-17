/**
 * This Kernel is implements a darcy flow problem using the full
 * Finite Element solution. It is included in this package for validation
 * purposes. It was slightly modified from the original DwarfElephantFEDarcy
 * Kernel in order to allow an easier usability within the Blackbox Wrapper of
 * OpenDA.
 */

///-------------------------------------------------------------------------
#ifndef DwarfElephantFEDarcyOPENDA_H
#define DwarfElephantFEDarcyOPENDA_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Kernel.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEDarcyOpenDA;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEDarcyOpenDA>();

///This Kernel is implements a darcy flow problem using the full Finite Element solution. It is included in this package for validation purposes. It was slightly modified from the original DwarfElephantFEDarcy Kernel in order to allow an easier usability within the Blackbox Wrapper of OpenDA.
class DwarfElephantFEDarcyOpenDA : public Kernel
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEDarcyOpenDA(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  std::vector<Real> _permeability;

  unsigned int _ID_first_block;
};

///-------------------------------------------------------------------------
#endif // DwarfElephantFEDarcyOpenDA_H
