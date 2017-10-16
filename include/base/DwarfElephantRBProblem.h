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

    virtual DwarfElephantRBAssembly & rbAssembly(unsigned int subdomain_id) { return *_rb_assembly[subdomain_id]; }

    virtual void newRBAssemblyArray(NonlinearSystemBase & nl);

    virtual MooseVariable & getVariable(THREAD_ID tid, const std::string & var_name) override;

//    std::vector<std::string> & getKernelNames(){return _kernel_names;}

  protected:
    NonlinearSystem * _nl_sys;

//    std::vector<std::string> _kernel_names;
    std::vector<DwarfElephantRBAssembly *> _rb_assembly;

};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBPROBLEM_H
