#ifndef DWARFELEPHANTFEFILEPOINTVALUES_H
#define DWARFELEPHANTFEFILEPOINTVALUES_H

#include "GeneralVectorPostprocessor.h"

// Forward Declarations
class DwarfElephantFEFilePointValues;

template <>
InputParameters validParams<DwarfElephantFEFilePointValues>();

///Extract the values of the physical points given by a text file.
class DwarfElephantFEFilePointValues : public GeneralVectorPostprocessor
{
public:
  DwarfElephantFEFilePointValues(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override {}

protected:
  const unsigned int _var_number;
  const System & _system;
  VectorPostprocessorValue & _values;
  std::string _file;
  std::vector<Real> _x_coordinates;
  std::vector<Real> _y_coordinates;
  std::vector<Real> _z_coordinates;
  unsigned int _num_points;

private:
  void fileParser();
};

#endif /* DWARFELEPHANTFEFILEPOINTVALUES_H */
