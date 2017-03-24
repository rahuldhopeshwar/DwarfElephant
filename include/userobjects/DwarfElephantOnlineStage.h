///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTONLINESTAGE_H
#define DWARFELEPHANTONLINESTAGE_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "GeneralUserObject.h"

///-------------------------------------------------------------------------
// Forward Declarations
class DwarfElephantOnlineStage;

template<>
InputParameters validParams<DwarfElephantOnlineStage>();

class DwarfElephantOnlineStage :
  public GeneralUserObject
{
  public:
    DwarfElephantOnlineStage(const InputParameters & params);

    void onlineStage();

    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTONLINESTAGE_H
