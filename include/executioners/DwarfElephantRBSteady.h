///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSTEADY_H
#define DWARFELEPHANTRBSTEADY_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Steady.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBSteady;

template<>
InputParameters validParams<DwarfElephantRBSteady>();

class DwarfElephantRBSteady :
  public Steady
{
  public:
    DwarfElephantRBSteady(const InputParameters & params);

    void execute() override;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSTEADY_H
