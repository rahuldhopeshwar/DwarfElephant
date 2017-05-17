///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBTRANSIENT_H
#define DWARFELEPHANTRBTRANSIENT_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Transient.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBTransient;

template<>
InputParameters validParams<DwarfElephantRBTransient>();

class DwarfElephantRBTransient :
  public Transient
{
  public:
    DwarfElephantRBTransient(const InputParameters & params);

    void execute() override;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBTRANSIENT_H
