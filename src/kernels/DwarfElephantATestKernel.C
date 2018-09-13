#include "DwarfElephantATestKernel.h"

//libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

template<>
InputParameters validParams<DwarfElephantATestKernel> ()
{
	InputParameters params  = validParams<Kernel>();
	params.addClassDescription("Implements an A matrix which makes use of a Empirically interpolated function");
        params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for  initializing the RB system");
	return params;
}

DwarfElephantATestKernel::DwarfElephantATestKernel(const InputParameters & parameters) :
    Kernel(parameters)
    
{
}

void DwarfElephantATestKernel::computeJacobian()
{		
        DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
        _local_ke.resize(ke.m(), ke.n());
        _local_ke.zero();        

	const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = getUserObject<DwarfElephantInitializeRBSystemSteadyState>("initial_rb_userobject");
        
        _local_ke.zero();

          for (_i = 0; _i < _test.size(); _i++)
	    for (_j = 0; _j < _phi.size(); _j++)  
        	for (_qp = 0; _qp < _qrule -> n_points(); _qp++)
		{
			_local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * _test[_i][_qp] * _phi[_j][_qp] * 1./sqrt(pow(_q_point[_qp](0) + 0.01,2) + pow(_q_point[_qp](1) + 0.01,2));
		}
	  ke += _local_ke;
	  if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() == 0)		
	  	_initialize_rb_system._fullFEnonAffineA -> add_matrix(_local_ke, _var.dofIndices());

}

Real DwarfElephantATestKernel::computeQpResidual()
{
	
	return 1.0;
}

Real DwarfElephantATestKernel::computeQpJacobian()
{
	return 0;
}
