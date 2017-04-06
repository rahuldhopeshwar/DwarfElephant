///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTOFFLINESTAGE_H
#define DWARFELEPHANTOFFLINESTAGE_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"

// MOOSE includes
#include "GeneralUserObject.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"
#include "NonlinearSystemBase.h"
#include "NodalBC.h"
#include "Assembly.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantRBClasses.h"
#include "DwarfElephantInitializeRBSystem.h"
#include "CacheStiffnessMatrix.h"


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
struct bnd_node_iterator;
class DwarfElephantOfflineStage;

template<>
InputParameters validParams<DwarfElephantOfflineStage>();

class DwarfElephantOfflineStage :
  public GeneralUserObject
{
  public:
    DwarfElephantOfflineStage(const InputParameters & params);

    void setInnerProductMatrix();
    void offlineStage();
    void setOnlineParameters();
    void transferAffineVectors();
    Real trainReducedBasis(const bool resize_rb_eval_data=true);
    Real truthSolve(int plot_solution);
    void truthAssembly();

    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;

  protected:
    bool _use_displaced;
    bool _store_basis_functions;
    bool _skip_matrix_assembly_in_rb_system;
    bool _skip_vector_assembly_in_rb_system;
    bool _compliant;
    bool _online_stage;

    std::string _system_name;
    std::string _residual_name;

    EquationSystems & _es;
    TransientNonlinearImplicitSystem & _sys;
    const DwarfElephantInitializeRBSystem & _initialize_rb_system;

    Function * _function;
    CacheStiffnessMatrix * _cache_stiffness_matrix;
    MooseMesh * _mesh_ptr;

    const std::set<SubdomainID> & _subdomain_ids;

    Real _mu_bar;
    unsigned int _online_N;
    std::vector<Real> _online_mu_parameters;

    RBParameters _rb_online_mu;

    MooseObjectWarehouse <NodalBC> _nodal_bcs;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTOFFLINESTAGE_H
