///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBEXECUTIONER_H
#define DWARFELEPHANTRBEXECUTIONER_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Steady.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBExecutioner;

template<>
InputParameters validParams<DwarfElephantRBExecutioner>();

class DwarfElephantRBExecutioner :
  public Steady
{
  public:
    DwarfElephantRBExecutioner(const InputParameters & params);

    void execute() override;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBEXECUTIONER_H
