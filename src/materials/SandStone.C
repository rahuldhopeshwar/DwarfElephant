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
 * permeability - http://www.geomore.com/porosity-and-permeability-2/
 *
 * viscosity - https://en.wikipedia.org/wiki/Brine
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
    _thermal_conductivity(declareProperty<Real>("thermal_conductivity")),
    _permeability(declareProperty<Real>("permeability")),
    _dynamic_viscosity(declareProperty<Real>("dynamic_viscosity")),
    _fluid_density(declareProperty<Real>("fluid_density")),
    _gravity(declareProperty<RealVectorValue>("gravity"))

{
}

///-------------------------------------------------------------------------
void
SandStone::computeQpProperties()
{
  _thermal_conductivity[_qp] = 2.5;                  // [W/(m K)]
  _permeability[_qp] = 1.e-6;                        // [m²]
  _dynamic_viscosity[_qp] =  0.001145;               // [kg/(m s)]
  _fluid_density[_qp] = 1;                        // [g/m³]
  _gravity[_qp] = RealVectorValue(0, -9.81, 0);       // [m/s²]
}
