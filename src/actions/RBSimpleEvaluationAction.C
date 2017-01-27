#include "RBSimpleEvaluationAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

template<>
InputParameters validParams<RBSimpleEvaluationAction>()
{
  InputParameters params = validParams<Action>();
  return params;
}

RBSimpleEvaluationAction::RBSimpleEvaluationAction(InputParameters params) :
    Action(params)
{
}

void
RBSimpleEvaluationAction::act()
{
}

