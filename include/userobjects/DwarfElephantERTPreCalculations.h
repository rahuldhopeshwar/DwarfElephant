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

///-------------------------------------------------------------------------
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

    DenseVector<Real> & computeGemeotricFactor();
    DenseVector<Real> & computeDistance(DenseVector<Real> & vec1, DenseVector<Real> & vec2);

//--------------------------------PROTECTED---------------------------------
  protected:

    /* Attributes */
    const std::vector<ExecFlagType> & _exec_flags;

    DenseVector<Real> _position_A_electrode;
    DenseVector<Real> _position_B_electrode;
    DenseVector<Real> _position_M_electrode;
    DenseVector<Real> _position_N_electrode;

    DenseVector<Real> _geometric_factor;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTERTPRECALCULATIONS_H
