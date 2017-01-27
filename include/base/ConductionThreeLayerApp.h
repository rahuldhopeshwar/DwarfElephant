/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#ifndef CONDUCTIONTHREELAYERAPP_H
#define CONDUCTIONTHREELAYERAPP_H

#include "MooseApp.h"

class ConductionThreeLayerApp;

template<>
InputParameters validParams<ConductionThreeLayerApp>();

class ConductionThreeLayerApp : public MooseApp
{
public:
  ConductionThreeLayerApp(InputParameters parameters);
  virtual ~ConductionThreeLayerApp();
  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax& syntax, ActionFactory & action_factory);
};

#endif /* CONDUCTIONTHREELAYERAPP_H */
