/* This UserObject is required to initialitze the RB system structure
 * and transfer for the steady state case.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTREDUCEDTOFULLSTATE_H
#define DWARFELEPHANTREDUCEDTOFULLSTATE_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/xdr_cxx.h"

// MOOSE includes
#include "GeneralUserObject.h"

///-------------------------------------------------------------------------
class DwarfElephantReducedToFullState;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<DwarfElephantReducedToFullState>();

///This UserObject is required to initialitze the RB system structure and transfer for the steady state case.
class DwarfElephantReducedToFullState :
  public GeneralUserObject
{

//----------------------------------PUBLIC----------------------------------
  public:
    DwarfElephantReducedToFullState(const InputParameters & params);

    /* Methods */

    // Initializes the RB System.
    virtual void initialize() override;

    // Method not used in this UserObject.
    virtual void execute() override;

    // Method not used in this UserObject.
    virtual void finalize() override;

    void read_reduced_state(const std::string & _file_base, const bool _read_binary_data);

    void check_file_exists(const std::string & _file_name);

//--------------------------------PROTECTED---------------------------------
  protected:

    /* Attributes */
    std::string _file_base;
    std::vector<Number> _state_vector_reduced;

    bool _read_binary_data;
};
///-------------------------------------------------------------------------
#endif // DWARFELEPHANTREDUCEDTOFULLSTATE_H
