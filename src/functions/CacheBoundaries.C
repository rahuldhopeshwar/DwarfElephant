/**
 * Function is used to cache the boundary conrtibutions for the stiffness
 * matrices and the load vectors.
 */

///---------------------------------INCLUDES--------------------------------
#include "CacheBoundaries.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<CacheBoundaries>()
{
  InputParameters params = validParams<Function>();
  return params;
}

///-------------------------------CONSTRUCTOR-------------------------------
CacheBoundaries::CacheBoundaries(const InputParameters & parameters) :
    Function(parameters)
{}

///-------------------------------------------------------------------------
Real
CacheBoundaries::value(Real /*t*/, const Point & /*p*/)
{
  return 0;
}

void
CacheBoundaries::cacheSubdomainResidual(numeric_index_type i, Real value, unsigned int subdomain)
{
  _cached_residual_subdomain_contribution_rows[subdomain].push_back(i);
  _cached_residual_subdomain_contribution_vals[subdomain].push_back(value);
}

void
CacheBoundaries::setCachedSubdomainResidual(NumericVector<Number> & _residual, unsigned int subdomain)
{
  _residual.close();

  for (unsigned int i = 0; i < _cached_residual_subdomain_contribution_vals[subdomain].size(); ++i)
    _residual.set(_cached_residual_subdomain_contribution_rows[subdomain][i], _cached_residual_subdomain_contribution_vals[subdomain][i]);
}

void
CacheBoundaries::cacheStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value)
{
  _cached_jacobian_contribution_rows.push_back(i);
  _cached_jacobian_contribution_cols.push_back(j);
  _cached_jacobian_contribution_vals.push_back(value);
}

void
CacheBoundaries::resizeSubdomainMatrixCaches(unsigned int subdomains)
{
  _cached_jacobian_subdomain_contribution_rows.resize(subdomains);
  _cached_jacobian_subdomain_contribution_cols.resize(subdomains);
  _cached_jacobian_subdomain_contribution_vals.resize(subdomains);
}

void
CacheBoundaries::resizeSubdomainVectorCaches(unsigned int subdomains)
{
  _cached_residual_subdomain_contribution_rows.resize(subdomains);
  _cached_residual_subdomain_contribution_vals.resize(subdomains);
}

void
CacheBoundaries::cacheSubdomainStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value, unsigned int subdomain)
{
  _cached_jacobian_subdomain_contribution_rows[subdomain].push_back(i);
  _cached_jacobian_subdomain_contribution_cols[subdomain].push_back(j);
  _cached_jacobian_subdomain_contribution_vals[subdomain].push_back(value);
}

void
CacheBoundaries::setCachedStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian)
{
  _jacobian.close();
  _jacobian.zero_rows(_cached_jacobian_contribution_rows);

  for (unsigned int i = 0; i < _cached_jacobian_contribution_vals.size(); ++i)
    _jacobian.set(_cached_jacobian_contribution_rows[i],
                  _cached_jacobian_contribution_cols[i],
                  _cached_jacobian_contribution_vals[i]);
}

void
CacheBoundaries::setCachedSubdomainStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian, unsigned int subdomain)
{
  _jacobian.close();
  _jacobian.zero_rows(_cached_jacobian_subdomain_contribution_rows[subdomain]);

  for (unsigned int i = 0; i < _cached_jacobian_subdomain_contribution_vals[subdomain].size(); ++i)
  {
    _jacobian.set(_cached_jacobian_subdomain_contribution_rows[subdomain][i],
                  _cached_jacobian_subdomain_contribution_cols[subdomain][i],
                  _cached_jacobian_subdomain_contribution_vals[subdomain][i]);
  }
}
