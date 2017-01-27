#ifndef DWARF_ELEPHANTAPP_H
#define DWARF_ELEPHANTAPP_H

#include "MooseApp.h"

class DwarfElephantApp;

template<>
InputParameters validParams<DwarfElephantApp>();

class DwarfElephantApp : public MooseApp
{
public:
  DwarfElephantApp(InputParameters parameters);
  virtual ~DwarfElephantApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* DWARF_ELEPHANTAPP_H */
