//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"

// MOOSE includes
#include "MortarConstraint.h"
#include "DisplacedProblem.h"
#include "NonlinearSystemBase.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"

// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
}

class NonlinearSystemBase;
class DisplacedProblem;

class DwarfElephantInitializeRBSystemSteadyState;
class DwarfElephantInitializeRBSystemTransient;
class DwarfElephantRBMortarConstraint;

template <>
InputParameters validParams<DwarfElephantRBMortarConstraint>();

class DwarfElephantRBMortarConstraint : public MortarConstraint
{
public:
  DwarfElephantRBMortarConstraint(const InputParameters & parameters);

  /* Methods */
   virtual void computeJacobian(Moose::MortarType mortar_type) override;
   virtual void computeResidual(Moose::MortarType mortar_type) override;
   virtual void initialSetup() override;
   // virtual void computeOutput();

  // Using declarations necessary to pull in computeResidual with different parameter list and avoid
  // hidden method warning
  //using MortarConstraintBase::computeResidual;

  // Using declarations necessary to pull in computeJacobian with different parameter list and avoid
  // hidden method warning
  //using MortarConstraintBase::computeJacobian;

protected:
  // virtual Real computeQpResidual(Moose::MortarType mortar_type) override;
  // virtual Real computeQpJacobian(Moose::ConstraintJacobianType jacobian_type,
  //                                unsigned int jvar) = 0;
  // virtual void computeQpOutput();

  /*Attributes*/
  FEProblemBase & _fe_problem;
  bool _use_displaced;
  bool _matrix_seperation_according_to_subdomains;
  bool _time_matrix_seperation_according_to_subdomains;
  bool _vector_seperation_according_to_subdomains;
  bool _compute_output;

  std::string _simulation_type;

  unsigned int _ID_first_block;
  unsigned int _ID_Aq;
  unsigned int _ID_Mq;
  unsigned int _ID_Fq;
  unsigned int _ID_Oq;

  Real _output_volume;

  DenseVector<Number> _local_out;

  EquationSystems & _es;

  const DwarfElephantInitializeRBSystemSteadyState * _initialize_rb_system;
  const DwarfElephantInitializeRBSystemTransient * _initialize_rb_system_transient;

};
