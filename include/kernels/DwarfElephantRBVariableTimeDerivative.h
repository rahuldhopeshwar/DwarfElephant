/**
 * This Kernel is implements a constant radiogenic heat production using the
 * Reduced Basis solution. It is required if you do want the radiogenic heat
 * production stay fixed for all parameters.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBVARIABLETIMEDERIVATIVE_H
#define DWARFELEPHANTRBVARIABLETIMEDERIVATIVE_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantRBTimeDerivative.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBVariableTimeDerivative;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBVariableTimeDerivative>();

///This Kernel is implements a constant radiogenic heat production using the  Reduced Basis solution. It is required if you do want the radiogenic heat production stay fixed for all parameters.
class DwarfElephantRBVariableTimeDerivative : public DwarfElephantRBTimeDerivative
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBVariableTimeDerivative(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBVARIABLETIMEDERIVATIVE_H
