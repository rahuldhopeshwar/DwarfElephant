#include "CacheStiffnessMatrix.h"

template<>
InputParameters validParams<CacheStiffnessMatrix>()
{
  InputParameters params = validParams<Function>();
  return params;
}

CacheStiffnessMatrix::CacheStiffnessMatrix(const InputParameters & parameters) :
    Function(parameters)
{}

Real
ExampleFunction::value(Real /*t*/, const Point & p)
{
  return 0;
}

void
CacheStiffnessMatrix::cacheStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value)
{
  _cached_jacobian_subdomain_contribution_rows.push_back(i,i);
  _cached_jacobian_subdomain_contribution_cols.push_back(j,j);
  _cached_jacobian_subdomain_contribution_vals.push_back(value);
}
