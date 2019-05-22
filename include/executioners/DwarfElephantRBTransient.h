#ifndef DWARFELEPHANTRBTRANSIENT_H
#define DWARFELEPHANTRBTRANSIENT_H

#include "Transient.h"

class DwarfElephantRBTransient;

template <>
InputParameters validParams<DwarfElephantRBTransient>();

class DwarfElephantRBTransient: public Transient
{
public:
  DwarfElephantRBTransient(const InputParameters & parameters);

  // virtual void init() override;

  virtual void execute() override;
  // virtual void init() override;
};

#endif // DWARFELEPHANTRBTRANSIENT_H
