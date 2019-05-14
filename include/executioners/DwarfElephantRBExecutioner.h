/**
 * This Executioner class is required to execute the RB method directly over
 * the libMesh system. It is for both the steady state and transient case.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBEXECUTIONER_H
#define DWARFELEPHANTRBEXECUTIONER_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "Steady.h"
#include "TimeStepper.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantRBExecutioner;

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBExecutioner>();

///This Executioner class is required to execute the RB method directly over the libMesh system. It is for both the steady state and transient case.
class DwarfElephantRBExecutioner :
  public Steady
{
//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantRBExecutioner(const InputParameters & params);

    /*Methods*/
    void execute() override;
    virtual Real computeDT() {return 0.0;};

//--------------------------------PROTECTED---------------------------------
  protected:
    /*Attributes*/
    std::string _simulation_type;
    bool _offline_stage;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBEXECUTIONER_H
