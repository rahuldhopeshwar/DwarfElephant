///-----------------------------------------------------------------------
#include "DwarfElephantRBTransientThetaExpansion.h"

DwarfElephantRBTransientThetaExpansion::DwarfElephantRBTransientThetaExpansion ()
  : Parent()
{}

Number
DwarfElephantRBTransientThetaExpansion::eval_IC_theta(unsigned int q,
                                                      const RBParameters & mu)
{
  if (q >= get_n_IC_terms())
    libmesh_error_msg("Error: We must have q < get_n_IC_terms in eval_IC_theta.");

  libmesh_assert(_IC_theta_vector[q]);

  return _IC_theta_vector[q]->evaluate( mu );
}

void
DwarfElephantRBTransientThetaExpansion::attach_IC_theta(RBTheta * theta_q_ic)
{
  libmesh_assert(theta_q_ic);

  _IC_theta_vector.push_back(theta_q_ic);
}
