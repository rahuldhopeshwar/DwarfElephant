#ifndef RBSIMPLEEVALUATIONACTION_H
#define RBSIMPLEEVALUATIONACTION_H

#include "Action.h"


//#include "RBConductionBase.h"


class RBSimpleEvaluationAction :
  public Action
{
public:
  RBSimpleEvaluationAction(InputParameters params);

  virtual void act() override;
};

template<>
InputParameters validParams<RBSimpleEvaluationAction>();

#endif //RBSIMPLEEVALUATIONACTION_H

