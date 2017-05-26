/**
 * Function is used to cache the boundary conrtibutions for the stiffness
 * matrices and the load vectors.
 */

///-------------------------------------------------------------------------
#ifndef CACHEBOUNDARIES_H
#define CACHEBOUNDARIES_H

///---------------------------------INCLUDES--------------------------------
// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"

// MOOSE includes
#include "Function.h"

///-------------------------------------------------------------------------
// Forward Declarations
namespace libMesh
{
  template <typename T> class SparseMatrix;
  template <typename T> class NumericVector;
}

class CacheBoundaries;

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<CacheBoundaries>();

///-------------------------------------------------------------------------
class CacheBoundaries : public Function
{
public:

//----------------------------------PUBLIC----------------------------------
  CacheBoundaries(const InputParameters & parameters);

  /* Methods */
  virtual Real value(Real t, const Point & p) override;
  void cacheStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value);
  void cacheSubdomainStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value, unsigned int subdomain);
  void cacheSubdomainResidual(numeric_index_type i, Real value, unsigned int subdomain);
  void cacheResidual(numeric_index_type i, Real value);

  void setCachedSubdomainResidual(NumericVector<Number> & _residual, unsigned int subdomain);
  void setCachedResidual(NumericVector<Number> & _residual);
  void setCachedStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian);
  void setCachedSubdomainStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian, unsigned int subdomain);

  void resizeSubdomainMatrixCaches(unsigned int subdomains);
  void resizeSubdomainVectorCaches(unsigned int subdomains);

//--------------------------------PROTECTED---------------------------------
protected:

  /* Attributes */
  std::vector <numeric_index_type> _cached_jacobian_contribution_rows;
  std::vector <numeric_index_type> _cached_jacobian_contribution_cols;
  std::vector <Real> _cached_jacobian_contribution_vals;

  std::vector<std::vector <numeric_index_type>> _cached_jacobian_subdomain_contribution_rows;
  std::vector<std::vector <numeric_index_type>> _cached_jacobian_subdomain_contribution_cols;
  std::vector<std::vector <Real>> _cached_jacobian_subdomain_contribution_vals;

  std::vector<std::vector <numeric_index_type>> _cached_residual_subdomain_contribution_rows;
  std::vector<std::vector <Real>> _cached_residual_subdomain_contribution_vals;

  std::vector <numeric_index_type> _cached_residual_contribution_rows;
  std::vector <Real> _cached_residual_contribution_vals;
};

#endif //CACHEBOUNDARIES_H
