#ifndef RBASSEMBLY_H
#define RBASSEMBLY_H

#include "Material.h"

//Forward Declarations
class RBAssembly;

template<>
InputParameters validParams<RBAssembly>();

class RBAssembly : public Material
{
public:
  RBAssembly(const InputParameters & parameters);

  virtual void computeProperties();
  void constructJacobianAndResidual();
};

#endif //RBASSEMBLY_H
