#ifndef DWARFELEPHANTRBELEMENTALVARIABLEVALUE_H
#define DWARFELEPHANTRBELEMENTALVARIABLEVALUE_H

#include "GeneralUserObject.h"

// Forward Declarations
class DwarfElephantRBElementalVariableValue;

class MooseMesh;

namespace libMesh
{
class Elem;
}


template <>
InputParameters validParams<DwarfElephantRBElementalVariableValue>();

class DwarfElephantRBElementalVariableValue : public GeneralUserObject
{
public:
  DwarfElephantRBElementalVariableValue(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}
  void assignElement(const std::vector<std::vector<NumericVector <Number> *> > _outputs);

protected:
  MooseMesh & _mesh;
  std::vector<unsigned int > _elementid;
  std::string _simulation_type;
};

#endif // DWARFELEPHANTRBELEMENTALVARIABLEVALUE_H
