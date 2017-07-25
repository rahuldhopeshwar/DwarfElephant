#ifndef DWARFELEPHANTRBASSEMBLY_H
#define DWARFELEPHANTRBASSEMBLY_H

// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"

//#include "Assembly.h"
#include "MooseTypes.h"

// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class NumericVector;
}

// MOOSE Forward Declares
class SystemBase;
//class Assembly;


class DwarfElephantRBAssembly //: public Assembly
{
public:
  DwarfElephantRBAssembly(SystemBase & sys, THREAD_ID tid);
  virtual ~DwarfElephantRBAssembly();

  void cacheStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value);
//  void cacheMassMatrixContribution(numeric_index_type i, numeric_index_type j, Real value);
  void cacheResidual(numeric_index_type i, Real value);

  void setCachedResidual(NumericVector<Number> & _residual);
  void setCachedStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian);
//  void setCachedMassMatrixContributions(SparseMatrix<Number> & _jacobian);

  void clearCachedStiffnessMatrixContributions();
//  void clearCachedMassMatrixContributions();


protected:
  SystemBase & _sys;
  THREAD_ID _tid;

  std::vector <numeric_index_type> _cached_jacobian_contribution_rows;
  std::vector <numeric_index_type> _cached_jacobian_contribution_cols;
  std::vector <Real> _cached_jacobian_contribution_vals;

//  std::vector <numeric_index_type> _cached_mass_contribution_rows;
//  std::vector <numeric_index_type> _cached_mass_contribution_cols;
//  std::vector <Real> _cached_mass_contribution_vals;

  std::vector <numeric_index_type> _cached_residual_contribution_rows;
  std::vector <Real> _cached_residual_contribution_vals;
};

#endif /* DWARFELEPHANTRBASSEMBLY_H */
