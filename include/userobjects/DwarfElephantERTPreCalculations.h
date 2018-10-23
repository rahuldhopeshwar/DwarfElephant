///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTERTPRECALCULATIONS_H
#define DWARFELEPHANTERTPRECALCULATIONS_H

///---------------------------------INCLUDES--------------------------------
//libMesh includes
#include "libmesh/equation_systems.h"

// MOOSE includes
#include "GeneralUserObject.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  class EquationSystems;
}

class DwarfElephantERTPreCalculations;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantERTPreCalculations>();

///Currently not used
class DwarfElephantERTPreCalculations :
  public GeneralUserObject
{

//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantERTPreCalculations(const InputParameters & params);

    /* Methods */
    // Initializes the RB System.
    virtual void initialize() override;

    // Method not used in this UserObject.
    virtual void execute() override;

    // Method not used in this UserObject.
    virtual void finalize() override;

    void setUpElectrodeVectors();
    NumericVector<Number> & computeGeometricFactor();
    NumericVector<Number> & computeDistance(NumericVector<Number> & vec1, NumericVector<Number> & vec2);

//--------------------------------PROTECTED---------------------------------
  protected:

    /* Attributes */
    const std::vector<ExecFlagType> & _exec_flags;

    unsigned int _n_electrodes;

    std::vector<Real> _position_A_electrode;
    std::vector<Real> _position_B_electrode;
    std::vector<Real> _position_M_electrode;
    std::vector<Real> _position_N_electrode;

    UniquePtr<NumericVector<Number>> _A_electrode;
    UniquePtr<NumericVector<Number>> _B_electrode;
    UniquePtr<NumericVector<Number>> _M_electrode;
    UniquePtr<NumericVector<Number>> _N_electrode;
    UniquePtr<NumericVector<Number>> _dist;
    UniquePtr<NumericVector<Number>> _geometric_factor;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTERTPRECALCULATIONS_H
