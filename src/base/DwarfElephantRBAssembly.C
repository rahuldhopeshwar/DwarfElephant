/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "DwarfElephantRBAssembly.h"

//// MOOSE includes
//#include "SubProblem.h"
//#include "ArbitraryQuadrature.h"
//#include "SystemBase.h"
//#include "MooseTypes.h"
//#include "MooseMesh.h"
//#include "MooseVariable.h"
//#include "MooseVariableScalar.h"
//#include "XFEMInterface.h"
//
//// libMesh
//#include "libmesh/coupling_matrix.h"
//#include "libmesh/dof_map.h"
//#include "libmesh/elem.h"
//#include "libmesh/equation_systems.h"
//#include "libmesh/fe_interface.h"
//#include "libmesh/node.h"
//#include "libmesh/quadrature_gauss.h"
//#include "libmesh/sparse_matrix.h"
//#include "libmesh/tensor_value.h"
//#include "libmesh/vector_value.h"

DwarfElephantRBAssembly::DwarfElephantRBAssembly(SystemBase & sys, THREAD_ID tid)
  : Assembly(sys, tid),
    _sys(sys),
    _tid(tid)
{
}

DwarfElephantRBAssembly::~DwarfElephantRBAssembly()
{
}
