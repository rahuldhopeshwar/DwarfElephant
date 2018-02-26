/**
 * This Output is required to couple the DwarfElephant package over a fork
 * interface to the OpenDA framework.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTOPENDAOUTPUT_H
#define DWARFELEPHANTOPENDAOUTPUT_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "FileOutput.h"

///-------------------------------------------------------------------------
// Forward declerations
class DwarfElephantOpenDAOutput;

///----------------------------INPUT PARAMETERS-----------------------------
template <>
InputParameters validParams<DwarfElephantOpenDAOutput>();

///-------------------------------------------------------------------------
class DwarfElephantOpenDAOutput : public FileOutput
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantOpenDAOutput(const InputParameters & parameters);

  /*Methods*/
  virtual void output(const ExecFlagType & type) override;

//--------------------------------PROTECTED---------------------------------
protected:
  /*Attributes*/
  std::string _system_name;
  VariableName _var_name;
  THREAD_ID _tid;
  PostprocessorName _postprocessor_name;
};

///-------------------------------------------------------------------------
#endif /* DWARFELEPHANTOPENDAOUTPUT_H */
