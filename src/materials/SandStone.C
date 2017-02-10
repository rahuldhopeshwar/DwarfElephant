/**
 * SandStone inherits from the MOOSE class Material and overrides
 * computeQpProperties.
 * Within the Class SandStone the typical rock properties of a sandstone are
 * stored.
 *
 * Attention: For a general use of the RB objects please use the string
 * "RBParameterNumber" for the parameters you want to include in the
 * RB Method. You can define RB Parameters for several properties and
 * store them in this file. But please make sure that all RB parameters
 * you do not want to include in the actual analysis are commented out.
 */

/* List of rock properties and corresponding references:
 *
 * thermal conductivity (lambda):
 *   http://www.minersoc.org/pages/Archive-CM/Volume_33/33-1-131.pdf
 *
 */

///---------------------------------INCLUDES--------------------------------
// MOOSE includes (DwarfElephant package)
#include "SandStone.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<SandStone>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Saves the rock properties of a typical \
                              sandstone");
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
SandStone::SandStone(const InputParameters & parameters) :
    Material(parameters),

    // Thermal conductivity
    _lambda(declareProperty<Real>("conductivity")),

    // RB parameters
    _rb_lambda(declareProperty<RealVectorValue>("mu0"))
{
}

///-------------------------------------------------------------------------
void
SandStone::computeQpProperties()
{
  // Thermal conductivity always in [W/(m K)]
  _lambda[_qp]= 2.5;

  // RB parameters
  _rb_lambda[_qp]={2.5, 3.0};
}
