#ifndef DWARFELEPHANTRBPOINTVALUE_H
#define DWARFELEPHANTRBPOINTVALUE_H

#include "GeneralUserObject.h"

// Forward Declarations
class DwarfElephantRBPointValue;

class MooseMesh;

namespace libMesh
{
class Elem;
}


template <>
InputParameters validParams<DwarfElephantRBPointValue>();

class DwarfElephantRBPointValue : public GeneralUserObject
{
public:
  DwarfElephantRBPointValue(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}
  void assignPoint(const std::vector<std::vector<NumericVector <Number> *> > _outputs);

protected:
  MooseMesh & _mesh;
  Point _point;
  unsigned int _outputid;
  std::string _simulation_type;
};

#endif // DWARFELEPHANTRBPOINTVALUE_H
