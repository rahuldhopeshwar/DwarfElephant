/* This UserObject is required to initialitze the RB system structure
 * and transfer for the 4D Variational Data Assimilation.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTINITIALIZERBSYSTEM4DVar_H
#define DWARFELEPHANTINITIALIZERBSYSTEM4DVar_H

///---------------------------------INCLUDES--------------------------------
// //libMesh includes
// #include "libmesh/equation_systems.h"
// #include "libmesh/sparse_matrix.h"
// #include "libmesh/petsc_matrix.h"
// #include "libmesh/petsc_vector.h"
//
// // MOOSE includes
// #include "GeneralUserObject.h"
// #include "DisplacedProblem.h"
// #include "MooseMesh.h"
//
// // MOOSE includes (DwarfElephant package)
// #include "DwarfElephantSystem.h"
#include "DwarfElephantRBClasses4DVar.h"
///-------------------------------------------------------------------------
// Forward Declarations
// namespace libMesh
// {
//   class EquationSystems;
// //  class RBConstruction;
//   template <typename T> class SparseMatrix;
//   template <typename T> class PetscMatrix;
//   template <typename T> class PetscVector;
// }

// class MooseMesh;
// class DwarfElephantSystem;
// class DwarfElephantRBConstructionTransient;
// class DwarfElephantRBConstructionSteadyState;
// class DwarfElephantInitializeRBSystemTransient;
class DwarfElephantInitializeRBSystem4DVar;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantInitializeRBSystem4DVar>();

///This UserObject is required to initialitze the RB system structure and transfer for the transient case.
class DwarfElephantInitializeRBSystem4DVar :
  public DwarfElephantInitializeRBSystemTransient
{

//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantInitializeRBSystem4DVar(const InputParameters & params);

    /* Methods */

    void processParameters();
    void initializeOfflineStage();
    void readObservationData(std::string _data_file);

    // Method not used in this UserObject.
    // virtual void execute() override;
    void execute();

    // Method not used in this UserObject.
    // virtual void finalize() override;
    void finalize();

    // std::vector<std::vector<NumericVector <Number> *> > getOutputs() const;

//--------------------------------PROTECTED---------------------------------
  protected:
    unsigned int _n_qois;
    std::vector<Real> _qoi_weights;
    std::string _data_file_name;

    // /*Friend Classes*/
    // friend class DwarfElephantRBKernel;
    // friend class DwarfElephantRBDiracKernel;
    // friend class DwarfElephantRBTimeKernel;
    // friend class DwarfElephantRBNodalBC;
    // friend class DwarfElephantRBIntegratedBC;
    // friend class DwarfElephantOfflineOnlineStageTransient;
    // friend class DwarfElephantRBExecutioner;
    // friend class DwarfElephantDakotaOutput;
    // friend class DwarfElephantRBInitialCondition;
    // friend class DwarfElephantRBProblem;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTINITIALIZERBSYSTEM4DVar_H
