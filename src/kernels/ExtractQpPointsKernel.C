#include "ExtractQpPointsKernel.h"

template<>
InputParameters validParams<ExtractQpPointsKernel>()
{
  InputParameters params = validParams<Diffusion>();
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
ExtractQpPointsKernel::ExtractQpPointsKernel(const InputParameters & parameters) :
  Diffusion(parameters)
{
}

///----------------------------------PDEs-----------------------------------
// Definition of the necessary PDE in the weak formulation
void
ExtractQpPointsKernel::initialSetup()
{
  // NonlinearSystemBase & _nl_sys = _fe_problem.getNonlinearSystemBase();
  _qp_file.open("QpPoints.txt", std::ios::app);
}

Real
ExtractQpPointsKernel::computeQpResidual()
{
  return Diffusion::computeQpResidual();
}

Real
ExtractQpPointsKernel::computeQpJacobian()
{
//  if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
//     _qp_file << _q_point(_qp) << std::endl;

  return Diffusion::computeQpJacobian();
}
