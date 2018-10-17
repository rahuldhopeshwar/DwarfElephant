#ifndef EXTRACTQPPOINTSKERNEL_H
#define EXTRACTQPPOINTSKERNEL_H

// MOOSE includes
#include "Diffusion.h"
#include "NonlinearSystemBase.h"


// Forward Declarations
///-------------------------------------------------------------------------
class ExtractQpPointsKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<ExtractQpPointsKernel>();


class ExtractQpPointsKernel : public Diffusion
{
//----------------------------------PUBLIC----------------------------------
public:
  ExtractQpPointsKernel(const InputParameters & parameters);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual void initialSetup() override;

  std::ofstream _qp_file;
};
#endif
