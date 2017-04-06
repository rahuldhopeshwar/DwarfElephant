#include "CacheStiffnessMatrix.h"

template<>
InputParameters validParams<CacheStiffnessMatrix>()
{
  InputParameters params = validParams<Function>();
  params.addParam<bool>("test_bool",true,"");
  return params;
}

CacheStiffnessMatrix::CacheStiffnessMatrix(const InputParameters & parameters) :
    Function(parameters)
{}

Real
CacheStiffnessMatrix::value(Real /*t*/, const Point & p)
{
  return 0;
}

void
CacheStiffnessMatrix::cacheSubdomainResidual(numeric_index_type i, Real value, unsigned int subdomain)
{
  _cached_residual_subdomain_contribution_rows[subdomain].push_back(i);
  _cached_residual_subdomain_contribution_vals[subdomain].push_back(value);
}

void
CacheStiffnessMatrix::setCachedSubdomainResidual(NumericVector<Number> & _residual, unsigned int subdomain)
{
  _residual.close();

  for (unsigned int i = 0; i < _cached_residual_subdomain_contribution_vals[subdomain].size(); ++i)
    _residual.set(_cached_residual_subdomain_contribution_rows[subdomain][i], _cached_residual_subdomain_contribution_vals[subdomain][i]*-1);

//  _residual.close();
}

void
CacheStiffnessMatrix::cacheStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value)
{
  _cached_jacobian_contribution_rows.push_back(i);
  _cached_jacobian_contribution_cols.push_back(j);
  _cached_jacobian_contribution_vals.push_back(value);
}

void
CacheStiffnessMatrix::resizeSubdomainMatrixCaches(unsigned int subdomains)
{
  _cached_jacobian_subdomain_contribution_rows.resize(subdomains);
  _cached_jacobian_subdomain_contribution_cols.resize(subdomains);
  _cached_jacobian_subdomain_contribution_vals.resize(subdomains);
}

void
CacheStiffnessMatrix::resizeSubdomainVectorCaches(unsigned int subdomains)
{
  _cached_residual_subdomain_contribution_rows.resize(subdomains);
  _cached_residual_subdomain_contribution_vals.resize(subdomains);
}

void
CacheStiffnessMatrix::cacheSubdomainStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value, unsigned int subdomain)
{
  _cached_jacobian_subdomain_contribution_rows[subdomain].push_back(i);
  _cached_jacobian_subdomain_contribution_cols[subdomain].push_back(j);
  _cached_jacobian_subdomain_contribution_vals[subdomain].push_back(value);
}

void
CacheStiffnessMatrix::setCachedStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian)
{
  _jacobian.close();
  _jacobian.zero_rows(_cached_jacobian_contribution_rows);

  for (unsigned int i = 0; i < _cached_jacobian_contribution_vals.size(); ++i)
    _jacobian.set(_cached_jacobian_contribution_rows[i],
                  _cached_jacobian_contribution_cols[i],
                  _cached_jacobian_contribution_vals[i]);
}

void
CacheStiffnessMatrix::setCachedSubdomainStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian, unsigned int subdomain)
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
