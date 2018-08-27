#ifndef DWARFELEPHANTRBFILEPOINTVALUES_H
#define DWARFELEPHANTRBFILEPOINTVALUES_H

#include "GeneralUserObject.h"

// Forward Declarations
class DwarfElephantRBFilePointValues;

class MooseMesh;

namespace libMesh
{
class Elem;
class FEInterface;
}


template <>
InputParameters validParams<DwarfElephantRBFilePointValues>();

class DwarfElephantRBFilePointValues : public GeneralUserObject
{
public:
  DwarfElephantRBFilePointValues(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}
  void assignPoint(const std::vector<std::vector<NumericVector <Number> *> > _outputs,
    unsigned int _outputid, Point _point ,bool _insistOnSuccess = false);

protected:
  MooseMesh & _mesh;
  std::string _file;
  const unsigned int _var_number;
  std::string _simulation_type;
  std::vector<Real> _x_coordinates;
  std::vector<Real> _y_coordinates;
  std::vector<Real> _z_coordinates;
  unsigned int _num_points;

private:
  void fileParser();
};

#endif // DWARFELEPHANTRBFILEPOINTVALUES_H
