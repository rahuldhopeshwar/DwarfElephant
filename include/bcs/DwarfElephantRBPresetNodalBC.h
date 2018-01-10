/**
 * This BC is required to use the RB method for the Preset BC. Since all
 * necessary operations are defined in the DwarfElephantRBNodalBC class this
 * class is the same as the PresetNodal BC provided by MOOSE except that it
 * inherits from the DwarfElephantRBNodalBC instead of the NodalBC class.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBPRESETNODALBC_H
#define DWARFELEPHANTRBPRESETNODALBC_H

//---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBDirichletBC.h"
#include "MooseVariable.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBPresetNodalBC;

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBPresetNodalBC>();

///-------------------------------------------------------------------------
class DwarfElephantRBPresetNodalBC : public DwarfElephantRBNodalBC
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBPresetNodalBC(const InputParameters & parameters);

  /*Methods*/
  void computeValue(NumericVector<Number> & current_solution);
//--------------------------------PROTECTED---------------------------------
protected:
  /*Methods*/
  virtual Real computeQpResidual() override;
  virtual Real computeQpValue();

};

#endif /* DWARFELEPHANTRBPRESETNODALBC_H */
