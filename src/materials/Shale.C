/**
 * Shale inherits from the MOOSE class  Material and overrides
 * computeQpProperties.
 * Within the class Shale the typical rock properties of a shale are stored.
 *
 * Attention: For a general use of the RB objects please use the string
 * "RBParameterNumber" for the parameters you want to include in the
 * RB Method. You can define RB Parameters for several properties and
 * store them in this file. But please make sure that all RB parameters
 * you do not want to include in the actual analysis are commented out.
 */

/* List of rock properties and corresponding references:
 *
 * thermal conductivuty (lambda) -
 *   http://www.minersoc.org/pages/Archive-CM/Volume_33/33-1-131.pdf
 *
 */

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "Shale.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<Shale>()
{
  InputParameters params = validParams<Material>();

   return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
Shale::Shale(const InputParameters & parameters) :
    Material(parameters),

    // Thermal conductivity
    _lambda(declareProperty<Real>("conductivity")),

    //RB Parameters
    _rb_lambda(declareProperty<RealVectorValue>("mu1"))
{
}

///-------------------------------------------------------------------------
void
Shale::computeQpProperties()
{
  // Thermal conductivity always in [W/(m K)]
  _lambda[_qp]= 1.05;

  // RB Parameters
  _rb_lambda[_qp]= {1.05, 1.55};
}
