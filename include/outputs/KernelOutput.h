#ifndef KERNELOUTPUT_H
#define KERNELOUTPUT_H

// MOOSE includes
#include "AdvancedOutput.h"
#include "Console.h"
#include "FileOutput.h"

class KernelOutput;

template<>
InputParameters validParams<KernelOutput>();

class KernelOutput :
  public AdvancedOutput<FileOutput>

{
public:

  KernelOutput(const InputParameters & parameters);

  virtual void output(const ExecFlagType & type);

protected:

   THREAD_ID _tid;
   AuxiliarySystem * _aux_sys_ptr;
   unsigned int _n_aux_var;
   MooseVariable & _var;
   bool _nodal;
   const VariableValue & _u;
   NumericVector< Number >  & _residual;
};

#endif // KERNELOUTPUT
