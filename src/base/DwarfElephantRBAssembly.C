#include "DwarfElephantRBAssembly.h"

DwarfElephantRBAssembly::DwarfElephantRBAssembly(SystemBase & sys, THREAD_ID tid)
  : //Assembly(sys, tid),
    _sys(sys),
    _tid(tid)
{
}

DwarfElephantRBAssembly::~DwarfElephantRBAssembly()
{
}

void
DwarfElephantRBAssembly::cacheSubdomainResidual(numeric_index_type i, Real value, unsigned int subdomain)
{
  _cached_residual_subdomain_contribution_rows[subdomain].push_back(i);
  _cached_residual_subdomain_contribution_vals[subdomain].push_back(value);
}

void
DwarfElephantRBAssembly::cacheResidual(numeric_index_type i, Real value)
{
  _cached_residual_contribution_rows.push_back(i);
  _cached_residual_contribution_vals.push_back(value);
}

void
DwarfElephantRBAssembly::setCachedSubdomainResidual(NumericVector<Number> & _residual, unsigned int subdomain)
{
  _residual.close();

  for (unsigned int i = 0; i < _cached_residual_subdomain_contribution_vals[subdomain].size(); ++i)
    _residual.set(_cached_residual_subdomain_contribution_rows[subdomain][i], _cached_residual_subdomain_contribution_vals[subdomain][i]);
}

void
DwarfElephantRBAssembly::setCachedResidual(NumericVector<Number> & _residual)
{
  _residual.close();

  for (unsigned int i = 0; i < _cached_residual_contribution_vals.size(); ++i)
    _residual.set(_cached_residual_contribution_rows[i], _cached_residual_contribution_vals[i]); // / 0.07778);
}

void
DwarfElephantRBAssembly::cacheStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value)
{
  _cached_jacobian_contribution_rows.push_back(i);
  _cached_jacobian_contribution_cols.push_back(j);
  _cached_jacobian_contribution_vals.push_back(value);
}

void
DwarfElephantRBAssembly::resizeSubdomainStiffnessMatrixCaches(unsigned int subdomains)
{
  _cached_jacobian_subdomain_contribution_rows.resize(subdomains);
  _cached_jacobian_subdomain_contribution_cols.resize(subdomains);
  _cached_jacobian_subdomain_contribution_vals.resize(subdomains);
}

void
DwarfElephantRBAssembly::resizeSubdomainMassMatrixCaches(unsigned int subdomains)
{
  _cached_mass_subdomain_contribution_rows.resize(subdomains);
  _cached_mass_subdomain_contribution_cols.resize(subdomains);
  _cached_mass_subdomain_contribution_vals.resize(subdomains);
}

void
DwarfElephantRBAssembly::resizeSubdomainVectorCaches(unsigned int subdomains)
{
  _cached_residual_subdomain_contribution_rows.resize(subdomains);
  _cached_residual_subdomain_contribution_vals.resize(subdomains);
}

void
DwarfElephantRBAssembly::cacheSubdomainStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value, unsigned int subdomain)
{
  _cached_jacobian_subdomain_contribution_rows[subdomain].push_back(i);
  _cached_jacobian_subdomain_contribution_cols[subdomain].push_back(j);
  _cached_jacobian_subdomain_contribution_vals[subdomain].push_back(value);
}

void
DwarfElephantRBAssembly::cacheSubdomainMassMatrixContribution(numeric_index_type i, numeric_index_type j, Real value, unsigned int subdomain)
{
  _cached_mass_subdomain_contribution_rows[subdomain].push_back(i);
  _cached_mass_subdomain_contribution_cols[subdomain].push_back(j);
  _cached_mass_subdomain_contribution_vals[subdomain].push_back(value);
}

void
DwarfElephantRBAssembly::setCachedStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian)
{
  _jacobian.close();
  _jacobian.zero_rows(_cached_jacobian_contribution_rows);

  for (unsigned int i = 0; i < _cached_jacobian_contribution_vals.size(); ++i)
    _jacobian.set(_cached_jacobian_contribution_rows[i],
                  _cached_jacobian_contribution_cols[i],
                  _cached_jacobian_contribution_vals[i]);
}

void
DwarfElephantRBAssembly::setCachedSubdomainStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian, unsigned int subdomain)
{
  _jacobian.close();
  _jacobian.zero_rows(_cached_jacobian_subdomain_contribution_rows[subdomain]);

  for (unsigned int i = 0; i < _cached_jacobian_subdomain_contribution_vals[subdomain].size(); ++i)
  {
    _jacobian.set(_cached_jacobian_subdomain_contribution_rows[subdomain][i],
                  _cached_jacobian_subdomain_contribution_cols[subdomain][i],
                  _cached_jacobian_subdomain_contribution_vals[subdomain][i]);
  }

  clearCachedSubdomainStiffnessMatrixContributions(subdomain);
}

void
DwarfElephantRBAssembly::setCachedSubdomainMassMatrixContributions(SparseMatrix<Number> & _mass, unsigned int subdomain)
{
  _mass.close();
  _mass.zero_rows(_cached_mass_subdomain_contribution_rows[subdomain]);

  for (unsigned int i = 0; i < _cached_mass_subdomain_contribution_vals[subdomain].size(); ++i)
  {
    _mass.set(_cached_mass_subdomain_contribution_rows[subdomain][i],
              _cached_mass_subdomain_contribution_cols[subdomain][i],
              _cached_mass_subdomain_contribution_vals[subdomain][i]);
  }
}


void
DwarfElephantRBAssembly::clearCachedSubdomainStiffnessMatrixContributions(unsigned int subdomains)
{
    unsigned int orig_size = _cached_jacobian_contribution_rows.size();

    _cached_jacobian_subdomain_contribution_rows[subdomains].clear();
    _cached_jacobian_subdomain_contribution_cols[subdomains].clear();
    _cached_jacobian_subdomain_contribution_vals[subdomains].clear();

    // It's possible (though massively unlikely) that clear() will
    // change the capacity of the vectors, so let's be paranoid and
    // explicitly reserve() the same amount of memory to avoid multiple
    // push_back() induced allocations.  We reserve 20% more than the
    // original size that was cached to account for variations in the
    // number of BCs assigned to each thread (for when the Jacobian
    // contributions are computed threaded).
    _cached_jacobian_subdomain_contribution_rows[subdomains].reserve(1.2 * orig_size);
    _cached_jacobian_subdomain_contribution_cols[subdomains].reserve(1.2 * orig_size);
    _cached_jacobian_subdomain_contribution_vals[subdomains].reserve(1.2 * orig_size);
}
