/**
 * This UserObject implements the Offline and Online stage of the RB method.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTOFFLINEONLINESTAGESTEADYSTATE_H
#define DWARFELEPHANTOFFLINEONLINESTAGESTEADYSTATE_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"

// MOOSE includes
#include "GeneralUserObject.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "NonlinearSystemBase.h"
#include "NodalBC.h"
#include "Assembly.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBClassesSteadyState.h"
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantRBProblem.h"


///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
  template <typename T> class SparseMatrix;
}

class MooseMesh;
class NonlinearSystemBase;
class Assembly;
class DwarfElephantOfflineOnlineStageSteadyState;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantOfflineOnlineStageSteadyState>();

///-------------------------------------------------------------------------
class DwarfElephantOfflineOnlineStageSteadyState :
  public GeneralUserObject
{

//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantOfflineOnlineStageSteadyState(const InputParameters & params);


    /* Methods */
    void setAffineMatrices();
    void offlineStage();
    void setOnlineParameters();
    void transferAffineVectors();
    std::string getFileName();

    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;

//--------------------------------PROTECTED---------------------------------
  protected:

    /* Attributes */
    bool _use_displaced;
    bool _store_basis_functions;
    bool _skip_matrix_assembly_in_rb_system;
    bool _skip_vector_assembly_in_rb_system;
    bool _compliant;
    bool _offline_stage;
    bool _online_stage;
    bool _offline_error_bound;
    bool _output_file;
    bool _compute_output;

    std::string _system_name;

    EquationSystems & _es;
    TransientNonlinearImplicitSystem & _sys;
    const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system;

    MooseMesh * _mesh_ptr;

    const std::set<SubdomainID> & _subdomain_ids;

    Real _mu_bar;
    unsigned int _online_N;
    std::vector<Real> _online_mu_parameters;

    RBParameters _rb_online_mu;

    DwarfElephantRBProblem * _rb_problem;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTOFFLINEONLINESTAGESTEADYSTATE_H
