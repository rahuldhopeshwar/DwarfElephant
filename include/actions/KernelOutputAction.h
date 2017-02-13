#ifndef KERNELOUTPUTACTION_H
#define KERNELOUTPUTACTION_H

#include "Action.h"

class KernelOutputAction : public Action
{
public:
  KernelOutputAction(InputParameters params);

  virtual void act() override;
};

template<>
InputParameters validParams<KernelOutputAction>();

#endif //KERNELOUTPUTACTION_H
