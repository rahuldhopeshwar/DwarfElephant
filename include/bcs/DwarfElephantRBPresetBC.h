/**
 * This BC is required to use the RB method for the Preset BC. Since all
 * necessary operations are defined in the DwarfElephantRBPresetNodalBC class
 * this class is the same as the Preset BC provided by MOOSE except that it
 * inherits from the DwarfElephantRBPresetNodalBC instead of the PresetNodalBC
 * class.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBPRESETBC_H
#define DWARFELEPHANTRBPRESETBC_H

//---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBPresetNodalBC.h"

//-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBPresetBC;

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBPresetBC>();

///This BC is required to use the RB method for the Preset BC. Since all necessary operations are defined in the DwarfElephantRBPresetNodalBC class this class is the same as the Preset BC provided by MOOSE except that it inherits from the DwarfElephantRBPresetNodalBC instead of the PresetNodalBC class.
class DwarfElephantRBPresetBC : public DwarfElephantRBPresetNodalBC
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBPresetBC(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:
  /*Methods*/
  virtual Real computeQpValue() override;

  /*Attributes*/
  const Real & _value;
};

#endif /* DWARFELEPHANTRBPRESETBC_H */
