///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTREDUCEDTOFULLSTATE_H
#define DWARFELEPHANTREDUCEDTOFULLSTATE_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/xdr_cxx.h"

// MOOSE includes
#include "GeneralUserObject.h"

///-------------------------------------------------------------------------
class DwarfElephantReducedToFullState;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantReducedToFullState>();

///This UserObject is required to transfer the reduced back into the full state.
class DwarfElephantReducedToFullState :
  public GeneralUserObject
{

//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantReducedToFullState(const InputParameters & params);

    /* Methods */

    // Initializes the RB System.
    virtual void initialize() override;

    // Method not used in this UserObject.
    virtual void execute() override;

    // Method not used in this UserObject.
    virtual void finalize() override;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTREDUCEDTOFULLSTATE_H
