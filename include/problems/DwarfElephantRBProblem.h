/**
 * This Problem class is required to use the RB system instead of the FE
 * system.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBPROBLEM_H
#define DWARFELEPHANTRBPROBLEM_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "FEProblemBase.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantSystem.h"
#include "DwarfElephantRBAssembly.h"


///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBProblem;
class NonlinearSystem;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBProblem>();

///-------------------------------------------------------------------------
class DwarfElephantRBProblem :
  public FEProblemBase
{
//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantRBProblem(const InputParameters & params);

    virtual ~DwarfElephantRBProblem();

    /* Methods */
    virtual void init () override;

    virtual void solve () override;

    virtual void setInputParametersFEProblem(InputParameters & parameters) override;

    NonlinearSystem & getNonlinearSystem() override { return *_nl_sys; }

    virtual DwarfElephantRBAssembly & rbAssembly(unsigned int subdomain_id) { return *_rb_assembly[subdomain_id]; }

    virtual void newRBAssemblyArray(NonlinearSystemBase & nl);

    // virtual MooseVariable & getVariable(THREAD_ID tid, const std::string & var_name) override;

//--------------------------------PROTECTED---------------------------------
  protected:
    /* Attributes */
    std::shared_ptr<NonlinearSystem> _nl_sys;
    // In case you are using a MOOSE version that is older than August 11 please replace the
    // line above with this line:
    //NonlinearSystem * _nl_sys;

    std::vector<DwarfElephantRBAssembly *> _rb_assembly;

};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBPROBLEM_H
