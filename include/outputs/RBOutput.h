#ifndef RBOUTPUT_H
#define RBOUTPUT_H

// MOOSE includes
#include "AdvancedOutput.h"
#include "Output.h"
#include "PetscOutput.h"
#include "Console.h"
#include "MooseMesh.h"
#include "FEProblem.h"
#include "FileMesh.h"
#include "MooseApp.h"
#include "NonlinearSystem.h"
#include "Assembly.h"
#include "GeneralizedPlaneStrainUserObject.h"
#include "ComputeJacobianBlocksThread.h"


// libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"

// RB classes
#include "libmesh/rb_construction.h"
#include "libmesh/rb_evaluation.h"

// RB classes simplified
//#include "RBConductionBase.h"
#include "RBStructuresP1ThetaEqualMu.h"
#include "rb_classes.h"
#include "RBKernel.h"
#include "RBSimpleConstruction.h"

// Forward declarations
namespace libMesh
{
  class EquationSystems;

  class RBConstruction;
  class RBEvaluation;

  class RBClasses;
  }

//class RBConductionBase;
class MooseMesh;
class JacobianBlock;
class RBOutput;

template<>
InputParameters validParams<RBOutput>();

  struct test : ElemAssembly, NonlinearSystem
{
    test(FEProblem & problem, const std::string & name_in):
      NonlinearSystem(problem, name_in)
    {}

  virtual void interior_assembly(FEMContext & c)
  {
    RBKernel * _rb_ptr;
//    c.get_elem_residual() = _rb_ptr-> &_local_ke;
  }

  friend class RBOutput;
  friend class RBClasses;
};


//struct test: NonlinearSystem
//  test(FEProblem & problem, const std::string & name_in):
//    NonlinearSystem(problem, name_in)
//    {}
//{
////  virtual void interior_assembly();
//////  {
////////    MooseApp * _app_ptr = getMooseApp();
////////    MooseObject * _moose_ptr;
////////    FEProblem * _problem_test_ptr; // = getParam<FEProblem *>("_fe_problem");
////////    NonlinearSystem * _non_sys_ptr;
////// //   _problem_test_ptr = getParam<FEProblem *>("_fe_problem");
////////    _non_sys_ptr = &_problem_test_ptr->getNonlinearSystem();
////////    NumericVector<Number> & _residual = _non_sys_ptr->residualVector(Moose::KT_ALL);
//////  }
//};

class RBOutput :
  public AdvancedOutput<FileOutput>

{
public:

  RBOutput(const InputParameters & parameters);

  virtual void output(const ExecFlagType & type);

  void constructJacobian();
  virtual void initRBSystem();
//  virtual void RBOffline();
//  virtual void RBOnline();

protected:

   std::string parameters_filename;

   bool _offline_stage;
   bool _online_stage;
   bool _store_basis_functions;

   unsigned int _online_N;
   Real _online_mu0_parameters;

  NonlinearSystem * _non_sys_ptr;
  AuxiliarySystem * _aux_sys_ptr;
  const KernelWarehouse & _kernel_warehouse;
  const MooseObjectWarehouse<KernelBase> & _no_time_kernels_ref;
  MooseSharedPointer<KernelBase> _kernel_base;

  MooseMesh * _mesh_ptr;
  MooseMesh * _mesh_ptr_test;


  const ExecuteMooseObjectWarehouse< UserObject > * _test_ptr;
  //KernelWarehouse _kernels;
  RBKernel * _rb_ptr;
  JacobianBlock * _jacobian_block_ptr;

  const THREAD_ID _tid;
  Assembly * _assembly_ptr;

  unsigned int _test_sys;
  DenseMatrix<Number> _jacobian;
  unsigned int _block_test;
  unsigned int _test;
  const CouplingMatrix * _ptr;

};

#endif // RBOUTPUT
