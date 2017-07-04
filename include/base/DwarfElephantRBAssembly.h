#ifndef DWARFELEPHANTRBASSEMBLY_H
#define DWARFELEPHANTRBASSEMBLY_H

#include "Assembly.h"

// MOOSE Forward Declares
class SystemBase;
class Assembly;


class DwarfElephantRBAssembly: public Assembly
{
public:
  DwarfElephantRBAssembly(SystemBase & sys, THREAD_ID tid);
  virtual ~DwarfElephantRBAssembly();

  void cacheStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value);
  void cacheSubdomainStiffnessMatrixContribution(numeric_index_type i, numeric_index_type j, Real value, unsigned int subdomain);
  void cacheSubdomainMassMatrixContribution(numeric_index_type i, numeric_index_type j, Real value, unsigned int subdomain);
  void cacheSubdomainResidual(numeric_index_type i, Real value, unsigned int subdomain);
  void cacheResidual(numeric_index_type i, Real value);

  void setCachedSubdomainResidual(NumericVector<Number> & _residual, unsigned int subdomain);
  void setCachedResidual(NumericVector<Number> & _residual);
  void setCachedStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian);
  void setCachedSubdomainStiffnessMatrixContributions(SparseMatrix<Number> & _jacobian, unsigned int subdomain);
  void setCachedSubdomainMassMatrixContributions(SparseMatrix<Number> & _mass, unsigned int subdomain);

  void resizeSubdomainStiffnessMatrixCaches(unsigned int subdomains);
  void resizeSubdomainMassMatrixCaches(unsigned int subdomains);
  void resizeSubdomainVectorCaches(unsigned int subdomains);

protected:
  SystemBase & _sys;
  THREAD_ID _tid;

  std::vector <numeric_index_type> _cached_jacobian_contribution_rows;
  std::vector <numeric_index_type> _cached_jacobian_contribution_cols;
  std::vector <Real> _cached_jacobian_contribution_vals;

  std::vector<std::vector <numeric_index_type>> _cached_jacobian_subdomain_contribution_rows;
  std::vector<std::vector <numeric_index_type>> _cached_jacobian_subdomain_contribution_cols;
  std::vector<std::vector <Real>> _cached_jacobian_subdomain_contribution_vals;

  std::vector<std::vector <numeric_index_type>> _cached_mass_subdomain_contribution_rows;
  std::vector<std::vector <numeric_index_type>> _cached_mass_subdomain_contribution_cols;
  std::vector<std::vector <Real>> _cached_mass_subdomain_contribution_vals;

  std::vector<std::vector <numeric_index_type>> _cached_residual_subdomain_contribution_rows;
  std::vector<std::vector <Real>> _cached_residual_subdomain_contribution_vals;

  std::vector <numeric_index_type> _cached_residual_contribution_rows;
  std::vector <Real> _cached_residual_contribution_vals;
};

#endif /* DWARFELEPHANTRBASSEMBLY_H */
