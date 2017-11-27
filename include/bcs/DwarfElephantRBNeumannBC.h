/**
 * This BC is required to use the RB method for the Neumann BC. Since all
 * necessary operations are defined in the DwarfElephantRBIntegratedBC class this
 * class is the same as the Neumann BC provided by MOOSE except that it
 * inherits from the DwarfElephantRBIntegratedBC instead of the IntegratedBC class.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPAHNTRBNEUMANNBC_H
#define DWARFELEPHANTRBNEUMANNBC_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBIntegratedBC.h"

///-------------------------------------------------------------------------
// Foward Declarations
class DwarfElephantRBNeumannBC;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBNeumannBC>();

///-------------------------------------------------------------------------
class DwarfElephantRBNeumannBC : public DwarfElephantRBIntegratedBC
{
//----------------------------------PUBLIC----------------------------------
public:

  DwarfElephantRBNeumannBC(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;

  /* Attributes */
  const Real & _value;
};

///-------------------------------------------------------------------------
#endif //DWARFELEPHANTRBNEUMANNBC_H
