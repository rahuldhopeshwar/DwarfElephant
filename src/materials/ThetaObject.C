/**
 * Theta is a RB material class. Within this class the theta object of
 * the differnt subdomains should be defined.
 */

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "ThetaObject.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<ThetaObject>()
{
  InputParameters params = validParams<Material>();
//  params.addClassDescription("Stores the theta object of the first subdomain.");
//  params.addRequiredParam<RBParameters>("mu","Mu value of the individual subdomain.");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
ThetaObject::ThetaObject(const InputParameters & parameters) :
    Material(parameters)
//    RBParameters(),
//    _mu(get_value("mu"))
{
}

///-------------------------------------------------------------------------
void
ThetaObject::computeQpProperties()
{
}
