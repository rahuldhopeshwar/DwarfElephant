///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTEXECUTIONER_H
#define DWARFELEPHANTEXECUTIONER_H

///---------------------------------INCLUDES--------------------------------

// MOOSE includes
#include "Steady.h"
#include "BlockRestrictable.h"


///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantExecutioner;

template<>
InputParameters validParams<DwarfElephantExecutioner>();

class DwarfElephantExecutioner :
  public Steady,
  public BlockRestrictable
{
  public:
    DwarfElephantExecutioner(const InputParameters & params);

    void execute() override;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTEXECUTIONER_H
