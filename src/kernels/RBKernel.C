/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBKernel inherits from the Kernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * stiffness matrix is needed separated in its subdomain contributions. In
 * addition it overwrites the function computeResidual.
 */

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/threads.h"
#include "libmesh/quadrature.h"

////MOOSE includes
#include "Assembly.h"
#include "MooseVariable.h"
#include "Problem.h"
#include "SubProblem.h"
#include "SystemBase.h"

//MOOSE includes (DwarfElephant package)
#include "RBKernel.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<RBKernel>()
{
  InputParameters params = validParams<Kernel>();

  params.addClassDescription("Overwrites the function computeJacobian. This is required because for the RB method the stiffness matrix needs to be saved in its subdomain contributions.");
  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addParam<unsigned int>("ID_first_block", 0, "ID of the first block in the mesh");
  params.addParam<unsigned int>("ID_Aq", 0, "ID of the current stiffness matrix");
  params.addParam<unsigned int>("ID_Mq", 0, "ID of the current mass matrix");
  params.addParam<unsigned int>("ID_Fq", 0, "ID of the current stiffness matrix");
  params.addParam<bool>("matrix_seperation_according_to_subdomains", true, "Tells whether the stiffness matrix is separated according to the subdomain_ids");
  params.addParam<bool>("time_matrix_seperation_according_to_subdomains", true, "Tells whether the mass matrix is separated according to the subdomain_ids");
  params.addParam<bool>("vector_seperation_according_to_subdomains", false, "Tells whether the load vector is separated according to the subdomain_ids");
//  params.addRequiredParam<Real>("max_x","Maximum extension of the volume of interest in x-direction.");
//  params.addRequiredParam<Real>("min_x","Minimum extension of the volume of interest in x-direction.");
//  params.addRequiredParam<Real>("max_y","Maximum extension of the volume of interest in y-direction.");
//  params.addRequiredParam<Real>("min_y","Minimum extension of the volume of interest in y-direction.");
//  params.addRequiredParam<Real>("max_z","Maximum extension of the volume of interest in z-direction.");
//  params.addRequiredParam<Real>("min_z","Minimum extension of the volume of interest in z-direction.");

  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
RBKernel::RBKernel(const InputParameters & parameters) :
    Kernel(parameters),
    _use_displaced(getParam<bool>("use_displaced")),
    _matrix_seperation_according_to_subdomains(getParam<bool>("matrix_seperation_according_to_subdomains")),
    _time_matrix_seperation_according_to_subdomains(getParam<bool>("time_matrix_seperation_according_to_subdomains")),
    _vector_seperation_according_to_subdomains(getParam<bool>("vector_seperation_according_to_subdomains")),
    _simulation_type(getParam<std::string>("simulation_type")),
    _ID_first_block(getParam<unsigned int>("ID_first_block")),
    _ID_Aq(getParam<unsigned int>("ID_Aq")),
    _ID_Mq(getParam<unsigned int>("ID_Mq")),
    _ID_Fq(getParam<unsigned int>("ID_Fq")),
//    _max_x(getParam<Real>("max_x")),
//    _min_x(getParam<Real>("min_x")),
//    _max_y(getParam<Real>("max_y")),
//    _min_y(getParam<Real>("min_y")),
//    _max_z(getParam<Real>("max_z")),
//    _min_z(getParam<Real>("min_z")),
    _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es())
//    _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject"))

{
}

///-------------------------------------------------------------------------
void
RBKernel::initialSetup()
{
//  if(_initialize_rb_system._exec_flags[0] != EXEC_INITIAL)
//    mooseError("The initialization of the RB system has to be executed on 'initial'. "
//               "You defined a wrong state in your 'execute_on' line in the Inputfile. "
//               "Please, correct your settings.");

  mooseInfo("For performing the reduced basis method a seperation of the stiffness matrix and the load vector according to "
            "the theta values is necessary. Therefore, the algorithm needs an ID for the matrices and vectors. The default "
            "setting seperates both the matrices and the vectors into the subdomain contributions. By performing the "
            "Kernels block wise and specify the ID in the inputfile any other seperation is also possible. For the seperation "
            "a default value of zero is assumed. In case your first volume ID is unequal to zero you have to define the "
            "'first_block_ID' in the input file. Due to the seperation the occurence of segmentation faults is likely. If "
            "a segementation fault occurs right after 'quiet mode?' check whether you: 1. defined the correct "
            "first_block_ID, 2. used the correct RBStructes header file in the RBClasses class.");

//  _output_volume = (_max_x - _min_x) * (_max_y - _min_y) * (_max_z - _min_z);
}

void
RBKernel::computeResidual()
{
  if(_vector_seperation_according_to_subdomains)
    _ID_Fq = _current_elem->subdomain_id() - _ID_first_block;
  
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

//  _local_out.resize(re.size());
//  _local_out.zero();

  Point _centroid = _current_elem->centroid();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); _i++)
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      {
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();
      //_console << "weights: " << _JxW[_qp] << std::endl;
      //_console << "residual: " << computeQpResidual() << std::endl;
      }

  //_console << *_fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0").rhs << std::endl;
  re += _local_re;


//  if ((_min_x <= _centroid(0)) && (_centroid(0) <= _max_x) &&
//      (_min_y <= _centroid(1)) && (_centroid(1) <= _max_y))
//    for (_i = 0; _i < _test.size(); _i++)
//      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//      _local_out(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual() / _output_volume;

  if(_simulation_type == "steady")  // SteadyState
  {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
      // Add the calculated vectors to the vectors from the RB system.
      //if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        _initialize_rb_system._residuals[_ID_Fq] -> add_vector(_local_re, _var.dofIndices());
  }
  else if (_simulation_type == "transient") // Transient
  {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
      // Add the calculated vectors to the vectors from the RB system.
      if (_fe_problem.getNonlinearSystemBase().computingInitialResidual())
        _initialize_rb_system._residuals[_ID_Fq] -> add_vector(_local_re, _var.dofIndices());
  }

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (const auto & var : _save_in)
      var->sys().solution().add_vector(_local_re, var->dofIndices());
  }
}

void
RBKernel::computeJacobian()
{
  if(_matrix_seperation_according_to_subdomains)
    _ID_Aq = _current_elem->subdomain_id() - _ID_first_block;

  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  precalculateJacobian();
  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _phi.size(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();



  ke += _local_ke;
  
  if(_simulation_type == "steady")  // Steady State
  {
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
    // Add the calculated matrices to the Aq matrices from the RB system.
    if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
        _initialize_rb_system._jacobian_subdomain[_ID_Aq] -> add_matrix(_local_ke, _var.dofIndices());
   }
     
  else if(_simulation_type == "transient") // Transient
  {
    const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
    if(_initialize_rb_system._offline_stage)
    // Add the calculated matrices to the Aq matrices from the RB system.
    if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)
    {
        _initialize_rb_system._jacobian_subdomain[_ID_Aq] -> add_matrix(_local_ke, _var.dofIndices());
        // Add the mass matrix to the RB System
        computeMassMatrix();
    }
  }
  
 if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i=0; i<rows; i++)
      diag(i) = _local_ke(i,i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (const auto & var : _diag_save_in)
      var->sys().solution().add_vector(diag, var->dofIndices());
  }
}

void
RBKernel::computeMassMatrix()
{
  if(_time_matrix_seperation_according_to_subdomains)
    _ID_Mq = _current_elem->subdomain_id() - _ID_first_block;

  DenseMatrix<Number> & me = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_me.resize(me.m(), me.n());
  _local_me.zero();

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _phi.size(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _local_me(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpMassMatrix();



  me += _local_me;
  const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemTransient>("initial_rb_userobject");
  _initialize_rb_system._mass_matrix_subdomain[_ID_Mq] -> add_matrix(_local_me, _var.dofIndices());

}
///----------------------------------PDEs-----------------------------------
// For the PDEs zero is implemented, since this Kernel shall be used for any
// RB problem. The problem specific PDEs are implemented in separate Kernels.

Real
RBKernel::computeQpJacobian()
{
  return 0;
}

Real
RBKernel::computeQpResidual()
{
  return 0;
}

Real
RBKernel::computeQpMassMatrix()
{
  return 0;
}
