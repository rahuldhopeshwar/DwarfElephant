///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBSYSTEM_H
#define DWARFELEPHANTRBSYSTEM_H

///---------------------------------INCLUDES--------------------------------
#include <iostream>
#include <fstream>

//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/getpot.h"

// MOOSE includes
#include "GeneralUserObject.h"
#include "DisplacedProblem.h"
#include "BlockRestrictable.h"
#include "MooseMesh.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBClasses.h"


///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
}

class MooseMesh;
class DwarfElephantRBSystem;

template<>
InputParameters validParams<DwarfElephantRBSystem>();

class DwarfElephantRBSystem :
  public GeneralUserObject,
  public BlockRestrictable

{
  public:
    DwarfElephantRBSystem(const InputParameters & params);

    void prepareDataStructuresRB();
    void offlineStage();
    void onlineStage();
    void performRBSystem();

    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;


  protected:
    bool _use_displaced;
    bool _skip_matrix_assembly;
    bool _skip_vector_assembly;
    bool _offline_stage;
    bool _online_stage;
    bool _store_basis_functions;

    unsigned int _online_N;

    Real _online_mu;

    NumericVector<Number> * _rhs_vector;
    NumericVector<Number> * _residual_non_time_vector;
    SparseMatrix<Number> * _system_matrix;

    std::string _system_name;
    std::string _parameters_filename;

    EquationSystems & _es;
    TransientNonlinearImplicitSystem & _sys;

    MooseMesh * _mesh_ptr;

    DwarfElephantRBConstruction & _rb_con;

    std::string _file_name;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTRBSYSTEM_H
