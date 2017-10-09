#ifndef EXTRACTQPPOINTS_H
#define EXTRACTQPPOINTS_H

// MOOSE includes
#include "Material.h"

// Forward Declarations
class ExtractQpPoints;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<ExtractQpPoints>();

///-------------------------------------------------------------------------
class ExtractQpPoints : public Material
{

//----------------------------------PUBLIC----------------------------------
public:
  ExtractQpPoints(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:
/* Methods */
  virtual void computeQpProperties() override;
};

///-------------------------------------------------------------------------
#endif //EXTRACTQPPOINTS_H
