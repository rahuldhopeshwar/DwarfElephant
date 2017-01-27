#ifndef RBSIMPLECONSTRUCTION_H
#define RBSIMPLECONSTRUCTION_H

/**
 * In this class simple subclasses of the classes RBEvaluation and
 * RBConstruction are introduced.
 *
 * RBSimpleEvaluation: requires only the definition of the lower coercivity
 * constant. The value is here specified for a Conduction problem.
 *
 * RBSimpleConstruction: In order to construct the RB System with the
 * RBSimpleEvaluation subclass the method build_rb_evaluation needs to be
 * overriden.
 */

#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_construction.h"
#include "RBKernel.h"

#include "RBStructuresP1ThetaEqualMu.h"
//#include "assembly.h"

#include "Action.h"

class RBSimpleConstruction;

template<>
InputParameters validParams<RBSimpleConstruction>();

class RBSimpleConstruction : public Action

{
public:
  RBSimpleConstruction(InputParameters & parameters);

//
//    SimpleRBConstruction (EquationSystems & es,
//                        const std::string & name_in,
//                        const unsigned int number_in)
//    : Parent(es, name_in, number_in),
//      dirichlet_bc(UniquePtr<DirichletBoundary>())

  virtual void act() override;

};


/////--------------------------RBSIMPLECONSTRUCTION--------------------------
//class RBSimpleConstruction : public RBConstruction
//{
//
////---------------------------------PUBLIC----------------------------------
//public:
//
//  // Constructor
//  RBSimpleConstruction (EquationSystems & _es,
//                              const std::string & _name_in,
//                              const unsigned int _number_in)
//    : Parent (_es, _name_in, _number_in)
//  {}
//
//  // Destructor
//  virtual ~RBSimpleConstruction () { }
//
//  // Type of the system
//  typedef RBSimpleConstruction _sys_type;
//
//  // Type of the parent
//  typedef RBConstruction Parent;
//
//  // Initialize data structure
//  virtual void init_data()
//  {
//
//    // Initialize the data from RBConstruction
//    Parent::init_data();
//
//    // Setting of the rb_assembly_expansion for the RBConstruction object.
//    // The theta expansion itself was set in the RBEvaluation object.
//    set_rb_assembly_expansion(_rb_assembly_expansion);
//
//    // Definition of the inner product matrix
//    set_inner_product_assembly(_rb_assembly_expansion.A0_assembly);
//  }
//
//
//  RBP1ThetaEqualMuAssemblyExpansion _rb_assembly_expansion;
//};


#endif // RBSIMPLECONSTRUCTION_H
