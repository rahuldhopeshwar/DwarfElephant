/**
 * This BC is required to use the RB method for the Dirichlet BC. Since all
 * necessary operations are defined in the DwarfElephantRBNodalBC class this
 * class is the same as the Dirichlet BC provided by MOOSE except that it
 * inherits from the DwarfElephantRBNodalBC instead of the NodalBC class.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBDIRICHLETBC_H
#define DWARFELEPHANTRBDIRICHLETBC_H

///---------------------------------INCLUDES--------------------------------
#include "DwarfElephantRBNodalBC.h"

///-------------------------------------------------------------------------
class DwarfElephantRBDirichletBC;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBDirichletBC>();

///-------------------------------------------------------------------------
class DwarfElephantRBDirichletBC : public DwarfElephantRBNodalBC
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBDirichletBC(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;

  /* Attributes */
  const Real & _value;
};

///-------------------------------------------------------------------------
#endif /* DWARFELEPHANTRBDIRICHLETBC_H */
