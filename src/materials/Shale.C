/// The material file Shale saves the rock properties of a typical shale.

/**
 * Attention: For a general use of the RB objects please use the string
 * "RBParameterNumber" for the parameters you want to include in the
 * RB Method. You can define RB Parameters for several properties and
 * store them in this file. But please make sure that all RB parameters
 * you do not want to include in the actual analysis are commented out.
 */

/* List of rock properties and corresponding references:
 *
 * thermal conductivuty (lambda) -
 * http://www.minersoc.org/pages/Archive-CM/Volume_33/33-1-131.pdf
 *
 */

//---------------------------------INCLUDE---------------------------------
#include "Shale.h"

//-------------------------------------------------------------------------
template<>

//----------------------------INPUT PARAMETERS-----------------------------
InputParameters validParams<Shale>()
{
  InputParameters params = validParams<Material>();

   return params;
}

//-----------------------------READ PARAMETERS-----------------------------
Shale::Shale(const InputParameters & parameters) :
    Material(parameters),

    /// Declare the material properties.  This returns references.

    /// Thermal conductivity
    _lambda(declareProperty<Real>("conductivity")),

    //RB Parameters
    _rb_lambda(declareProperty<RealVectorValue>("RBParameter0")),
    _online_lambda(declareProperty<Real>("OnlineConductivity"))
{
}

//-------------------------------------------------------------------------
void
Shale::computeQpProperties()
{
  /// Thermal conductivity always in [W/(m K)]
  _lambda[_qp]= 1.05;

  // RB Parameters
  _rb_lambda[_qp]= {1.05, 1.55};
  _online_lambda[_qp]=1.05;
}
