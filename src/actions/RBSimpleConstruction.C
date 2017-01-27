/**
 * In this class simple subclasses of the classes RBEvaluation and
 * RBConstruction are introduced.
 *
 * RBSimpleEvaluation: requires only the definition of the lower coercivity
 * constant. The value is here specified for a Conduction problem.
 *
 * RBSimpleConstruction: In order to construct the RB System with the
 * RBSimpleEvaluation subclass the method build_rb_evaluation needs to be
 * overriden.
 */

#include "RBSimpleConstruction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

template<>
InputParameters validParams<RBSimpleConstruction>()
{
  InputParameters params = validParams<Action>();
  return params;
}

RBSimpleConstruction::RBSimpleConstruction(InputParameters & parameters) :
    Action(parameters)
{
}
void
RBSimpleConstruction::act()
{
}
