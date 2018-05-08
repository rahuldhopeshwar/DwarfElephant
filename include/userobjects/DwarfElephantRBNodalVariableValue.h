#ifndef DWARFELEPHANTRBNODALVARIABLEVALUE_H
#define DWARFELEPHANTRBNODALVARIABLEVALUE_H

#include "GeneralUserObject.h"

// Forward Declarations
class DwarfElephantRBNodalVariableValue;

class MooseMesh;

namespace libMesh
{
class Node;
}


template <>
InputParameters validParams<DwarfElephantRBNodalVariableValue>();

class DwarfElephantRBNodalVariableValue : public GeneralUserObject
{
public:
  DwarfElephantRBNodalVariableValue(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}
  void assignNode(const std::vector<std::vector<NumericVector <Number> *> > _outputs);

protected:
  MooseMesh & _mesh;
  std::vector<unsigned int > _nodeid;
  const Real _scale_factor;
  std::string _simulation_type;
};

#endif // DWARFELEPHANTRBNODALVARIABLEVALUE_H
