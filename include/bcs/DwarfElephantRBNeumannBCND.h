/**
 * This BC is required to use the RB method for the dimensionless Neumann BC.
 * Since all necessary operations are defined in the DwarfElephantRBIntegratedBC class
 * this class is the same as the Neumann BC provided by MOOSE except that it
 * inherits from the DwarfElephantRBIntegratedBC instead of the IntegratedBC class.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPAHNTRBNEUMANNBCND_H
#define DWARFELEPHANTRBNEUMANNBCND_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBIntegratedBC.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBNeumannBCND;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBNeumannBCND>();

///-------------------------------------------------------------------------
class DwarfElephantRBNeumannBCND : public DwarfElephantRBIntegratedBC
{
//----------------------------------PUBLIC----------------------------------
public:

  DwarfElephantRBNeumannBCND(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;

  /* Attributes */
  const Real & _value;
  Real _u_ref;
  Real _l_ref;
};

///-------------------------------------------------------------------------
#endif //DWARFELEPHANTRBNEUMANNBC_H
