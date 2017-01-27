 /**
 * The structures are defined for an elliptic PDE with the following restrictions:
 *  1. The parametere dimension p is equal to one (an extension to p = 2 will
 *    be found in RBKernelSteadyStateP2ThetatEqualMu and an extension to
 *    p = 3 will be implemented in RBKernelSteadyStateP3ThetaEqualMu).
 * 2. Theta is equal to mu (for implementing other relationships,please
 *    follow the structure of these implementation for a general usability).
 *
 * The structures defined are:
 * 1. Theta --> parameter-dependent part of the PDE
 * 2. Aq --> stiffness matrix (parameter-independent)
 * 3. Fq --> load vector (parameter-independent)
 * 4. Output
 * 5. RBThetaExpansion
 * 6. RBAssemblyExpansion
 */

