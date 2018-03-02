#ifndef DWARFELEPHANTRBDIRACKERNEL_H
#define DWARFELEPHANTRBDIRACKERNEL_H

// Moose Includes
#include "DiracKernel.h"

// MOOSE includes (DwarfElephant package)
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"

// Forward Declarations
class DwarfElephantRBDiracKernel;

template <>
InputParameters validParams<DwarfElephantRBDiracKernel>();

class DwarfElephantRBDiracKernel : public DiracKernel
{
public:
  DwarfElephantRBDiracKernel(const InputParameters & parameters);

  virtual void computeResidual() override;
  virtual void computeJacobian() override;

protected:
  bool _matrix_seperation_according_to_subdomains;
  bool _vector_seperation_according_to_subdomains;

  std::string _simulation_type;

  unsigned int _ID_first_block;
  unsigned int _ID_Aq;
  unsigned int _ID_Fq;

  virtual Real computeQpResidual() override;
};

#endif // DWARFELEPHANTRBDIRACKERNEL_H
