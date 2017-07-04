#ifndef DWARFELEPHANTRBASSEMBLY_H
#define DWARFELEPHANTRBASSEMBLY_H

//#include "MooseArray.h"
#include "MooseTypes.h"
//#include "Assembly.h"

//// libMesh
//#include "libmesh/dense_matrix.h"
//#include "libmesh/dense_vector.h"
//#include "libmesh/enum_quadrature_type.h"
//#include "libmesh/fe_type.h"

//// libMesh forward declarations
//namespace libMesh
//{
//class DofMap;
//class CouplingMatrix;
//class Elem;
//template <typename T>
//class FEGenericBase;
//typedef FEGenericBase<Real> FEBase;
//class Node;
//template <typename T>
//class NumericVector;
//template <typename T>
//class SparseMatrix;
//
//typedef VectorValue<Real> RealVectorValue;
//typedef RealVectorValue RealGradient;
//
//typedef TensorValue<Real> RealTensorValue;
//typedef RealTensorValue RealTensor;
//}

// MOOSE Forward Declares
//class MooseMesh;
//class ArbitraryQuadrature;
class SystemBase;
//class MooseVariable;
//class XFEMInterface;
//typedef MooseArray<std::vector<Real>> VariablePhiValue;
//typedef MooseArray<std::vector<RealGradient>> VariablePhiGradient;
//typedef MooseArray<std::vector<RealTensor>> VariablePhiSecond;
class Assembly;


class DwarfElephantRBAssembly: public Assembly
{
public:
  DwarfElephantRBAssembly(SystemBase & sys, THREAD_ID tid);
  virtual ~DwarfElephantRBAssembly();
};

#endif /* DWARFELEPHANTRBASSEMBLY_H */
