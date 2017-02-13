#ifndef RBOUTPUT_H
#define RBOUTPUT_H

// libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/getpot.h"

// libMesh includes (RB package)
#include "libmesh/rb_construction.h"
#include "libmesh/rb_evaluation.h"


// MOOSE includes
#include "AdvancedOutput.h"
#include "Console.h"
#include "MooseMesh.h"
#include "FileMesh.h"

// MOOSE includes (DwarfElephant package)
#include "RBSimpleClasses.h"

// Forward declarations
namespace libMesh
{
  class EquationSystems;

  class RBConstruction;
  class RBEvaluation;

  class RBClasses;
  }

class MooseMesh;
class RBOutput;

template<>
InputParameters validParams<RBOutput>();

class RBOutput :
  public AdvancedOutput<FileOutput>

{
public:

  RBOutput(const InputParameters & parameters);

  virtual void output(const ExecFlagType & type);

  virtual void performRBSystem();
//  virtual void RBOffline();
//  virtual void RBOnline();

protected:

   std::string parameters_filename;

   bool _offline_stage;
   bool _online_stage;
   bool _store_basis_functions;

   unsigned int _online_N;
   Real _online_mu0_parameters;

   MooseMesh * _mesh_ptr;
//   THREAD_ID _tid;
//   AuxiliarySystem * _aux_sys_ptr;
//   unsigned int _n_aux_var;
//   MooseVariable & _var;
};

#endif // RBOUTPUT
