/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "Conduction.h"
template<>

InputParameters validParams<Conduction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  // in the case that lambda should be defined in the input file
  // params.addRequiredParam<Real>("lambda","thermal conductivity");
  return params;
}

Conduction::Conduction(const InputParameters & parameters) :
  Kernel(parameters),
  // in the case that lambda should be defined in the input file
  // _lambda(getParam<Real>("lambda"))

  // gets the thermal conductivity directly from the corresponding material file
  _lambda(getMaterialProperty<Real>("conductivity"))
{
}


// Definition of the necessary PDE in the weak formulation
Real
Conduction::computeQpResidual()
{
  return _lambda[_qp]*_grad_u[_qp] * _grad_test[_i][_qp];
}

Real
Conduction::computeQpJacobian()
{
  return _lambda[_qp]*_grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
