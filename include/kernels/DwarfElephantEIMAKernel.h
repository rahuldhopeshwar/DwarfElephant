#ifndef DWARFELEPHANTEIMAKERNEL_H
#define DWARFELEPHANTEIMAKERNEL_H

#include "Kernel.h"
#include "DisplacedProblem.h"
#include "NonlinearSystemBase.h"

#include "DwarfElephantInitializeRBSystemSteadyState.h"

#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"

//class DwarfElephantInitializeRBSystemSteadyState;

class DwarfElephantEIMAKernel;

template<>
InputParameters validParams<DwarfElephantEIMAKernel>();

class DwarfElephantEIMAKernel : public Kernel
{
	public:
	DwarfElephantEIMAKernel(const InputParameters & parameters);
        virtual void computeJacobian() override;
	
	protected:
	virtual Real computeQpResidual() override;
	virtual Real computeQpJacobian() override;
	
	std::vector<Number> _eim_values;
        int elem_number = 0;
};

#endif //DWARFELEPHANTEIMAKERNEL_H
