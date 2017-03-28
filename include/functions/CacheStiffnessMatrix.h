#ifndef CACHESTIFFNESSMATRIX_H
#define CACHESTIFFNESSMATRIX_H

#include "Function.h"

class CacheStiffnessMatrix;

template<>
InputParameters validParams<CacheStiffnessMatrix>();

class CacheStiffnessMatrix : public Function
{
public:
  CacheStiffnessMatrix(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) override;

protected:
  std::vector <numeric_index_type> _cached_jacobian_subdomain_contribution_rows;
  std::vector <numeric_index_type> _cached_jacobian_subdomain_contribution_cols;
  std::vector <Real> _cached_jacobian_subdomain_contribution_vals;
};

#endif //CACHESTIFFNESSMATRIX_H
