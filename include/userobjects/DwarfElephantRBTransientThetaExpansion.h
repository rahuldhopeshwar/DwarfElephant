///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTRBTRANSIENTTHETAEXPANSION_H
#define DWARFELEPHANTRBTRANSIENTTHETAEXPANSION_H

///---------------------------------INCLUDES--------------------------------
#include "libmesh/transient_rb_theta_expansion.h"
///-----------------------DWARFELEPHANTRBCONSTRUCTION-----------------------
class DwarfElephantRBTransientThetaExpansion : public TransientRBThetaExpansion
{

//---------------------------------PUBLIC-----------------------------------
public:

  // Constructor
  DwarfElephantRBTransientThetaExpansion ();

  typedef TransientRBThetaExpansion Parent;

  /**
   * Evaluate theta at the current parameter. Override
   * if the theta functions need to be treated differently
   * in subclasses.
   */
  virtual Number eval_IC_theta(unsigned int q, const RBParameters & mu);

  /**
   * Get Q_ic, the number of terms in the affine
   * expansion for the initial conditions.
   */
  virtual unsigned int get_n_IC_terms()
    { return cast_int<unsigned int>(_IC_theta_vector.size());}

  /**
   * Attach a pointer to a functor object that defines one
   * of the theta_q_ic terms.
   */
  virtual void attach_IC_theta(RBTheta * theta_q_ic);

private:
  /**
   * Vector storing the pointers to the RBTheta functors.
   */
  std::vector<RBTheta *> _IC_theta_vector;

};

#endif // DWARFELEPHANTRBTRANSIENTTHETAEXPANSION_H
