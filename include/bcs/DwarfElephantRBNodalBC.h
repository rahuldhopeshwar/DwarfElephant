/**
 * This BC is required to use the RB method as it is provided by the
 * RB libMesh package. The RBNodalBC inherits from the NodalBC class. It
 * overwrites the function computeJacobian because for the RB method the
 * stiffness matrix is needed separated in its subdomain contributions. In
 * addition it overwrites the function computeResidual.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBNODALBC_H
#define DWARFELEPHANTRBNODALBC_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/equation_systems.h"

// MOOSE includes
#include "NodalBC.h"
#include "MooseMesh.h"
#include "NonlinearSystemBase.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"
#include "DwarfElephantRBProblem.h"

///-------------------------------------------------------------------------
// Forward declarations
namespace libMesh
{
  class EquationSystems;
}

class MooseMesh;
class NonlinearSystemBase;

class DwarfElephantInitializeRBSystemSteadyState;
class DwarfElephantInitializeRBSystemTransient;
class DwarfElephantRBNodalBC;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantRBNodalBC>();

///-------------------------------------------------------------------------
class DwarfElephantRBNodalBC :
  public NodalBC
{

//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBNodalBC(const InputParameters & parameters);

  /* Methods */
  virtual void computeResidual() override;
  // for older Moose versions that still need the residual vector as an input
  // virtual void computeResidual(NumericVector<Number> & residual) override;
  virtual void computeJacobian() override;
  virtual void initialSetup() override;

//--------------------------------PROTECTED---------------------------------
protected:

  /* Methods */
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /* Attributes */
  bool _matrix_seperation_according_to_subdomains;
  bool _compute_output;

  std::string _simulation_type;

  std::vector<unsigned int> _ID_Fq;
  unsigned int _ID_Aq;
  unsigned int _ID_Mq;
  DwarfElephantRBProblem * _rb_problem;

  const DwarfElephantInitializeRBSystemSteadyState * _initialize_rb_system;
  const DwarfElephantInitializeRBSystemTransient * _initialize_rb_system_transient;
};

#endif /* DWARFELEPHANTRBNODALBC_H */
