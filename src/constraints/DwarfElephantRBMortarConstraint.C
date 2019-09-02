//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// libMesh includes
#include "libmesh/threads.h"
#include "libmesh/quadrature.h"

//MOOSE includes
#include "Assembly.h"
#include "MooseVariable.h"
#include "MooseVariableScalar.h"
#include "Problem.h"
#include "SubProblem.h"
#include "SystemBase.h"

//MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBMortarConstraint.h"

template <>
InputParameters
validParams<DwarfElephantRBMortarConstraint>()
{
  InputParameters params = validParams<MortarConstraint>();

  params.addClassDescription("RB object for mortar constraint");
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addParam<unsigned int>("ID_Aq", 0, "ID of the current stiffness matrix");
  params.addParam<unsigned int>("ID_Mq", 0, "ID of the current mass matrix");
  params.addParam<unsigned int>("ID_Fq", 0, "ID of the current load vector");
  params.addParam<unsigned int>("ID_Oq", 0, "ID of the current output vector");
  params.addParam<bool>("matrix_seperation_according_to_subdomains", true, "Tells whether the stiffness matrix is separated according to the subdomain_ids");
  params.addParam<bool>("time_matrix_seperation_according_to_subdomains", true, "Tells whether the mass matrix is separated according to the subdomain_ids");
  params.addParam<bool>("vector_seperation_according_to_subdomains", false, "Tells whether the load vector is separated according to the subdomain_ids");
  params.addParam<bool>("compute_output",false,"Determines whether an output function is used or not");

  return params;
}

DwarfElephantRBMortarConstraint::DwarfElephantRBMortarConstraint(const InputParameters & parameters)
  : MortarConstraint(parameters),
  _fe_problem(*getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
  _use_displaced(getParam<bool>("use_displaced")),
  _matrix_seperation_according_to_subdomains(getParam<bool>("matrix_seperation_according_to_subdomains")),
  _time_matrix_seperation_according_to_subdomains(getParam<bool>("time_matrix_seperation_according_to_subdomains")),
  _vector_seperation_according_to_subdomains(getParam<bool>("vector_seperation_according_to_subdomains")),
  _compute_output(getParam<bool>("compute_output")),
  _simulation_type(getParam<std::string>("simulation_type")),
  _ID_first_block(*_fe_problem.mesh().meshSubdomains().begin()),
  _ID_Aq(getParam<unsigned int>("ID_Aq")),
  _ID_Mq(getParam<unsigned int>("ID_Mq")),
  _ID_Fq(getParam<unsigned int>("ID_Fq")),
  _ID_Oq(getParam<unsigned int>("ID_Oq")),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es())
{

}

void
DwarfElephantRBMortarConstraint::initialSetup()
{
  if(_simulation_type == "steady")  // SteadyState
  {
    _initialize_rb_system = &getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");

    if(_initialize_rb_system->_exec_flags[0] != EXEC_INITIAL)
    mooseError("The UserObject 'DwarfElephantInitializeRBSystemSteadyState' has to be executed on 'initial'. "
               "You defined a wrong state in your 'execute_on' line in the input file. "
               "Please, correct your settings.");
  }
  else
  {
    _initialize_rb_system_transient = &getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");

    if(_initialize_rb_system_transient->_exec_flags[0] != EXEC_INITIAL)
    mooseError("The UserObject 'DwarfElephantInitializeRBSystemTransient' has to be executed on 'initial'. "
               "You defined a wrong state in your 'execute_on' line in the input file. "
               "Please, correct your settings.");
  }
}

void
DwarfElephantRBMortarConstraint::computeResidual(Moose::MortarType mortar_type)
{
   // if(_vector_seperation_according_to_subdomains)
   //   _ID_Fq = _current_elem->subdomain_id() - _ID_first_block;

  DenseVector<Number> & re = _assembly.residualBlock(_slave_var.number());
  unsigned int test_space_size = 0;

  switch (mortar_type)
  {
    case Moose::MortarType::Slave:
      prepareVectorTag(_assembly, _slave_var.number());
      //re = _assembly.residualBlock(_slave_var.number());
      test_space_size = _test_slave.size();
      _local_re.resize(re.size());
      _local_re.zero();
      break;

    case Moose::MortarType::Master:
      prepareVectorTagNeighbor(_assembly, _master_var.number());
      //re = _assembly.residualBlockNeighbor(_master_var.number());
      test_space_size = _test_master.size();
      _local_re.resize(re.size());
      _local_re.zero();
      break;

    case Moose::MortarType::Lower:
      mooseAssert(_var, "LM variable is null");
      prepareVectorTagLower(_assembly, _var->number());
      //re = _assembly.residualBlockLower(_var->number());
      test_space_size = _test.size();
      _local_re.resize(re.size());
      _local_re.zero();
      break;
  }

  for (_qp = 0; _qp < _qrule_msm->n_points(); _qp++)
    for (_i = 0; _i < test_space_size; _i++)
      _local_re(_i) += _JxW_msm[_qp] * _coord[_qp] * computeQpResidual(mortar_type);

  re += _local_re;

  // if(_simulation_type == "steady")  // SteadyState
  // {
  //   if (_ID_Fq >= _initialize_rb_system->_qf)
  //     mooseError("The number of load vectors you defined here is not matching the number of load vectors you specified in the RBClasses Class.");
  //
  //   if(_initialize_rb_system->_offline_stage)
  //     // Add the calculated vectors to the vectors from the RB system.
  //     if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
  //     {
  //       switch (mortar_type) {
  //         case Moose::MortarType::Slave:
  //           _initialize_rb_system->_residuals[_ID_Fq] -> add_vector(_local_re, _slave_var.dofIndices());
  //         break;
  //         case Moose::MortarType::Master:
  //           _initialize_rb_system->_residuals[_ID_Fq] -> add_vector(_local_re, _master_var.dofIndicesNeighbor());
  //         break;
  //         case Moose::MortarType::Lower:
  //           _initialize_rb_system->_residuals[_ID_Fq] -> add_vector(_local_re, _var->dofIndicesLower());
  //         break;
  //       }
  //     }
  // }
  // else if (_simulation_type == "transient") // Transient
  // {
  //   if (_ID_Fq >= _initialize_rb_system_transient->_qf)
  //     mooseError("The number of load vectors you defined here is not matching the number of load vectors you specified in the RBClasses Class.");
  //
  //   if(_initialize_rb_system_transient->_offline_stage)
  //     // Add the calculated vectors to the vectors from the RB system.
  //     if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
  //     {
  //       switch (mortar_type) {
  //         case Moose::MortarType::Slave:
  //           _initialize_rb_system_transient->_residuals[_ID_Fq] -> add_vector(_local_re, _slave_var.dofIndices());
  //         break;
  //         case Moose::MortarType::Master:
  //           _initialize_rb_system_transient->_residuals[_ID_Fq] -> add_vector(_local_re, _master_var.dofIndicesNeighbor());
  //         break;
  //         case Moose::MortarType::Lower:
  //           _initialize_rb_system_transient->_residuals[_ID_Fq] -> add_vector(_local_re, _var->dofIndicesLower());
  //         break;
  //       }
  //     }
  // }
}

void
DwarfElephantRBMortarConstraint::computeJacobian(Moose::MortarType mortar_type)
{
  size_t test_space_size = 0;
  typedef Moose::ConstraintJacobianType JType;
  typedef Moose::MortarType MType;
  std::array<JType, 3> jacobian_types;

  switch (mortar_type)
  {
    case MType::Slave:
      test_space_size = _slave_var.dofIndices().size();
      jacobian_types = {{JType::SlaveSlave, JType::SlaveMaster, JType::SlaveLower}};
      break;

    case MType::Master:
      test_space_size = _master_var.dofIndicesNeighbor().size();
      jacobian_types = {{JType::MasterSlave, JType::MasterMaster, JType::MasterLower}};
      break;

    case MType::Lower:
      test_space_size = _var ? _var->dofIndicesLower().size() : 0;
      jacobian_types = {{JType::LowerSlave, JType::LowerMaster, JType::LowerLower}};
      break;
  }

  auto & ce = _assembly.couplingEntries();
  for (const auto & it : ce)
  {
    MooseVariableFEBase & ivariable = *(it.first);
    MooseVariableFEBase & jvariable = *(it.second);

    unsigned int ivar = ivariable.number();
    unsigned int jvar = jvariable.number();

    switch (mortar_type)
    {
      case MType::Slave:
        if (ivar != _slave_var.number())
          continue;
        break;

      case MType::Master:
        if (ivar != _master_var.number())
          continue;
        break;

      case MType::Lower:
        if (!_var || _var->number() != ivar)
          continue;
        break;
    }

    std::array<size_t, 3> shape_space_sizes{{jvariable.dofIndices().size(),
                                             jvariable.dofIndicesNeighbor().size(),
                                             jvariable.dofIndicesLower().size()}};
    std::array<const VariablePhiValue *, 3> phis;
    std::array<const VariablePhiGradient *, 3> grad_phis;
    std::array<const VectorVariablePhiValue *, 3> vector_phis;
    std::array<const VectorVariablePhiGradient *, 3> vector_grad_phis;
    if (jvariable.isVector())
    {
      const auto & temp_var = static_cast<MooseVariableFE<RealVectorValue> &>(jvariable);
      vector_phis = {{&temp_var.phiFace(), &temp_var.phiFaceNeighbor(), &temp_var.phiLower()}};
      vector_grad_phis = {
          {&temp_var.gradPhiFace(), &temp_var.gradPhiFaceNeighbor(), &temp_var.gradPhiLower()}};
    }
    else
    {
      const auto & temp_var = static_cast<MooseVariableFE<Real> &>(jvariable);
      phis = {{&temp_var.phiFace(), &temp_var.phiFaceNeighbor(), &temp_var.phiLower()}};
      grad_phis = {
          {&temp_var.gradPhiFace(), &temp_var.gradPhiFaceNeighbor(), &temp_var.gradPhiLower()}};
    }

    for (MooseIndex(3) type_index = 0; type_index < 3; ++type_index)
    {
      prepareMatrixTagLower(_assembly, ivar, jvar, jacobian_types[type_index]);

      /// Set the proper phis
      if (jvariable.isVector())
      {
        _vector_phi = vector_phis[type_index];
        _vector_grad_phi = vector_grad_phis[type_index];
      }
      else
      {
        _phi = phis[type_index];
        _grad_phi = grad_phis[type_index];
      }

      for (_i = 0; _i < test_space_size; _i++)
        for (_j = 0; _j < shape_space_sizes[type_index]; _j++)
          for (_qp = 0; _qp < _qrule_msm->n_points(); _qp++)
            _local_ke(_i, _j) +=
                _JxW_msm[_qp] * _coord[_qp] * computeQpJacobian(jacobian_types[type_index], jvar);
      accumulateTaggedLocalMatrix();
    }
  }
}
