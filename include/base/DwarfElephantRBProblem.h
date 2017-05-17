///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBPROBLEM_H
#define DWARFELEPHANTRBPROBLEM_H

///---------------------------------INCLUDES--------------------------------

// MOOSE includes
#include "FEProblemBase.h"
#include "DwarfElephantSystem.h"


///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBProblem;
class NonlinearSystem;

template<>
InputParameters validParams<DwarfElephantRBProblem>();

class DwarfElephantRBProblem :
  public FEProblemBase
{
  public:
    DwarfElephantRBProblem(const InputParameters & params);

    virtual ~DwarfElephantRBProblem();

    virtual void solve () override;

    virtual void setInputParametersFEProblem(InputParameters & parameters) override;

    NonlinearSystem & getNonlinearSystem() override { return *_nl_sys; }

  protected:
    NonlinearSystem * _nl_sys;

};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBPROBLEM_H
