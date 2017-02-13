 /**
  * The structures are defined for an elliptic PDE with the following restrictions:
  *  1. The parameter dimension p is equal to one. (P1)
  *  2. Theta has three different values. (_3)
  *  3. Theta is equal to mu (for implementing other relationships,please
  *     follow the structure of these implementation for a general usability).
  *     (ThetaEqualMu)
  *  4. The problem is compliant. (Compliant)
  *
  * The structures defined are:
  * 1. Theta --> parameter-dependent part of the PDE
  * 2. Aq --> stiffness matrix (parameter-independent)
  * 3. Fq --> load vector (parameter-independent)
  * 4. Output
  * 5. RBThetaExpansion
  * 6. RBAssemblyExpansion
  */

///-------------------------------------------------------------------------

