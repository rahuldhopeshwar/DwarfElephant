/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBKernel inherits from the Kernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * whole stiffness matrix is needed and not only the diagonal entries.
 */

///-------------------------------------------------------------------------
#ifndef RBKERNEL_H
#define RBKERNEL_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"

// MOOSE includes
#include "Kernel.h"
#include "DisplacedProblem.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBClasses.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
}

class DwarfElephantRBConstruction;
class DisplacedProblem;
class RBKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<RBKernel>();

///-------------------------------------------------------------------------
class RBKernel : public Kernel
{

//----------------------------------PUBLIC----------------------------------
public:
  RBKernel(const InputParameters & parameters);

 /* Methods */
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void timestepSetup() override;

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  std::string _system_name;

  bool _use_displaced;

  EquationSystems & _es;
  TransientNonlinearImplicitSystem & _sys;
  DwarfElephantRBConstruction * _rb_con;

  std::vector <SparseMatrix<Number> * > _Aq_qa;
  SparseMatrix<Number> * _jacobian;
  };

///-------------------------------------------------------------------------
#endif //RBKERNEL_H
