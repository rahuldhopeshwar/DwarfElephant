///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBPROBLEM_H
#define DWARFELEPHANTRBPROBLEM_H

///---------------------------------INCLUDES--------------------------------

// MOOSE includes
#include "FEProblem.h"
#include "BlockRestrictable.h"


///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBProblem;

template<>
InputParameters validParams<DwarfElephantRBProblem>();

class DwarfElephantRBProblem :
  public FEProblem
{
  public:
    DwarfElephantRBProblem(const InputParameters & params);
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBPROBLEM_H
