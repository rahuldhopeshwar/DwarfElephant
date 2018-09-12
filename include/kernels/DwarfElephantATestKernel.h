#ifndef DWARFELEPHANTATESTKERNEL_H
#define DWARFELEPHANTATESTKERNEL_H

#include "Kernel.h"
#include "DisplacedProblem.h"
#include "NonlinearSystemBase.h"

#include "DwarfElephantInitializeRBSystemSteadyState.h"

#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"

//class DwarfElephantInitializeRBSystemSteadyState;

class DwarfElephantATestKernel;

template<>
InputParameters validParams<DwarfElephantATestKernel>();

class DwarfElephantATestKernel : public Kernel
{
	public:
	DwarfElephantATestKernel(const InputParameters & parameters);
        virtual void computeJacobian() override;
	
	protected:
	virtual Real computeQpResidual() override;
	virtual Real computeQpJacobian() override;
	
};

#endif //DWARFELEPHANTEIMAKERNEL_H
