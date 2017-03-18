/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBKernel inherits from the Kernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * stiffness matrix is needed separated in its subdomain contributions.
 */

///-------------------------------------------------------------------------
#ifndef RBKERNEL_H
#define RBKERNEL_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"

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
  template <typename T> class PetscMatrix;
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
  virtual void computeJacobian() override;
  virtual void initialSetup() override;
  virtual void timestepSetup() override;

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
  
  /*Attributes*/
  bool _use_displaced;

  EquationSystems & _es;
  DwarfElephantRBConstruction * _rb_con_ptr;
  
  const std::set<SubdomainID> & _block_ids;

  SparseMatrix<Number> * _jacobian_subdomain;
  };

///-------------------------------------------------------------------------
#endif //RBKERNEL_H
