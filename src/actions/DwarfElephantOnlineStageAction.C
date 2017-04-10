#include "DwarfElephantOnlineStageAction.h"

template<>
InputParameters validParams<DwarfElephantOnlineStageAction>()
{
  InputParameters params = validParams<Action>();
  return params;
}

DwarfElephantOnlineStageAction::DwarfElephantOnlineStageAction(InputParameters params) :
    Action(params)
{
}

void
DwarfElephantOnlineStageAction::act()
{
}

