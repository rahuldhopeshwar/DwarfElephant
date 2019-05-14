/**
 * This Kernel is implements a constant radiogenic heat production using the full
 * Finite Element solution. It is included in this package for validation
 * purposes.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTFECONSTANTSPECIFICSTORAGE_H
#define DWARFELEPHANTFECONSTANTSPECIFICSTORAGE_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "TimeDerivative.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantFEConstantSpecificStorage;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantFEConstantSpecificStorage>();

///This Kernel is implements a constant radiogenic heat production using the full Finite Element solution. It is included in this package for validation purposes.
class DwarfElephantFEConstantSpecificStorage : public TimeDerivative
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantFEConstantSpecificStorage(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  Real _specific_storage;
  Real _norm_value;
};

///-------------------------------------------------------------------------
#endif // DWARFELEPHANTFECONSTANTSPECIFICSTORAGE_H
