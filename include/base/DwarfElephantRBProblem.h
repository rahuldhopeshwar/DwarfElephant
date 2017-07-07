///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBPROBLEM_H
#define DWARFELEPHANTRBPROBLEM_H

///---------------------------------INCLUDES--------------------------------

// MOOSE includes
#include "FEProblemBase.h"
#include "DwarfElephantSystem.h"
#include "DwarfElephantRBAssembly.h"


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

    virtual DwarfElephantRBAssembly & rbAssembly(THREAD_ID tid) { return *_rb_assembly[tid]; }

    virtual void newRBAssemblyArray(NonlinearSystemBase & nl);
    //virtual void newAssemblyArray(NonlinearSystemBase & nl) override;

  protected:
    NonlinearSystem * _nl_sys;

    std::vector<DwarfElephantRBAssembly *> _rb_assembly;

};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBPROBLEM_H
