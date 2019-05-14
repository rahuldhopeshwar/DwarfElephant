/**
 * This Kernel is implements a constant radiogenic heat production using the
 * Reduced Basis solution. It is required if you do want the radiogenic heat
 * production stay fixed for all parameters.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBCONSTANTSPECIFICSTORAGE_H
#define DWARFELEPHANTRBCONSTANTSPECIFICSTORAGE_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "DwarfElephantRBTimeDerivative.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBConstantSpecificStorage;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBConstantSpecificStorage>();

///This Kernel is implements a constant radiogenic heat production using the  Reduced Basis solution. It is required if you do want the radiogenic heat production stay fixed for all parameters.
class DwarfElephantRBConstantSpecificStorage : public DwarfElephantRBTimeDerivative
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBConstantSpecificStorage(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  Real _specific_storage;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBCONSTANTSPECIFICSTORAGE_H
