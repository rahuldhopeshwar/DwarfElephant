/**
 * This Assembly class is required to cache the matrix and vector entries
 * introdued by the Dirichlet BCs with respect to the desired separation.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBASSEMBLY_H
#define DWARFELEPHANTRBASSEMBLY_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"

// MOOSE includes
#include "MooseTypes.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class NumericVector;
}

class SystemBase;

///-------------------------------------------------------------------------
class DwarfElephantRBAssembly
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantRBAssembly(SystemBase & sys, THREAD_ID tid);
  virtual ~DwarfElephantRBAssembly();

  /* Methods */
  void cacheJacobianContribution(numeric_index_type i, numeric_index_type j, Real value);
  void cacheMassMatrixContribution(numeric_index_type i, numeric_index_type j, Real value);
  void cacheResidual(numeric_index_type i, Real value);
  void cacheOutput(numeric_index_type i, Real value);

  void setCachedResidual(NumericVector<Number> & _residual);
  void setCachedOutput(NumericVector<Number> & _output);
  void setCachedJacobianContributions(SparseMatrix<Number> & _jacobian, bool _set_diagonal_entries);
  void setCachedMassMatrixContributions(SparseMatrix<Number> & _mass);

  void clearCachedJacobianContributions();
  void clearCachedMassMatrixContributions();
  void clearCachedResidualContributions();
  void clearCachedOutputContributions();

//--------------------------------PROTECTED---------------------------------
protected:
  /* Attributes */
  SystemBase & _sys;
  THREAD_ID _tid;

  std::vector <numeric_index_type> _cached_jacobian_contribution_rows;
  std::vector <numeric_index_type> _cached_jacobian_contribution_cols;
  std::vector <Real> _cached_jacobian_contribution_vals;

  std::vector <numeric_index_type> _cached_mass_contribution_rows;
  std::vector <numeric_index_type> _cached_mass_contribution_cols;
  std::vector <Real> _cached_mass_contribution_vals;

  std::vector <numeric_index_type> _cached_residual_contribution_rows;
  std::vector <Real> _cached_residual_contribution_vals;

  std::vector <numeric_index_type> _cached_output_contribution_rows;
  std::vector <Real> _cached_output_contribution_vals;
};

#endif /* DWARFELEPHANTRBASSEMBLY_H */
