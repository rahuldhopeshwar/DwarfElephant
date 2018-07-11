#ifndef DWARFELEPHANTRBINITIALCONDITION_H
#define DWARFELEPHANTRBINITIALCONDITION_H

#include "InitialCondition.h"

#include "DwarfElephantInitializeRBSystemTransient.h"

#include <string>

// Forward Declarations
class DwarfElephantInitializeRBSystemTransient;
class DwarfElephantRBInitialCondition;

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<DwarfElephantRBInitialCondition>();

class DwarfElephantRBInitialCondition : public InitialCondition
{
public:
  DwarfElephantRBInitialCondition(const InputParameters & parameters);

  virtual void initialSetup() override;

  virtual void compute() override;

protected:
  const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system;

  unsigned int _ID_first_block;
  unsigned int _ID_IC_q;
  unsigned int _ID_IC_q_for_shared_nodes;

  bool _vector_separation_according_to_subdomains;
  bool _shared_node_separation_according_to_subdomains;
};

#endif // DWARFELEPHANTRBINITIALCONDITION_H
