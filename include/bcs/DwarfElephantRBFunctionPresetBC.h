#ifndef DWARFELEPHANTRBFUNCTIONPRESETBC_H
#define DWARFELEPHANTRBFUNCTIONPRESETBC_H

#include "DwarfElephantRBPresetNodalBC.h"

//Forward Declarations
class DwarfElephantRBFunctionPresetBC;
class Function;

template<>
InputParameters validParams<DwarfElephantRBFunctionPresetBC>();

class DwarfElephantRBFunctionPresetBC : public DwarfElephantRBPresetNodalBC
{
public:
  DwarfElephantRBFunctionPresetBC(const InputParameters & parameters);

protected:

  virtual Real computeQpValue() override;

  Function & _func;
};

#endif //FUNCTIONPRESETBC_H
