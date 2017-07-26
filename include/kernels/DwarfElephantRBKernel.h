/**
 * This Kernel is required to use the RB method as it is provided by the
 * RB libMesh package. The RBKernel inherits from the Kernel class. It
 * overwrites the function computeJacobian because for the RB method the
 * stiffness matrix is needed separated in its subdomain contributions. In
 * addition it overwrites the function computeResidual.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBKERNEL_H
#define DWARFELEPHANTRBKERNEL_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"

// MOOSE includes
#include "Kernel.h"
#include "DisplacedProblem.h"
#include "NonlinearSystemBase.h"

// MOOSE includes (DwarfElephant package)
//#include "DwarfElephantRBClasses.h"
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
}

class NonlinearSystemBase;
class DisplacedProblem;

class DwarfElephantInitializeRBSystemSteadyState;
class DwarfElephantRBKernel;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBKernel>();

///-------------------------------------------------------------------------
class DwarfElephantRBKernel : public Kernel
{

//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBKernel(const InputParameters & parameters);

 /* Methods */
  virtual void computeJacobian() override;
  virtual void computeResidual() override;
  virtual void initialSetup() override;

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  /*Attributes*/
  bool _use_displaced;
  bool _matrix_seperation_according_to_subdomains;
  bool _time_matrix_seperation_according_to_subdomains;
  bool _vector_seperation_according_to_subdomains;

  std::string _simulation_type;

  unsigned int _ID_first_block;
  unsigned int _ID_Aq;
  unsigned int _ID_Mq;
  unsigned int _ID_Fq;

//  Real _max_x;
//  Real _min_x;
//  Real _max_y;
//  Real _min_y;
//  Real _max_z;
//  Real _min_z;
//  Real _output_volume;

//  DenseVector<Number> _local_out;

  EquationSystems & _es;
};

///-------------------------------------------------------------------------
#endif //DWARFELEPHANTRBKERNEL_H
