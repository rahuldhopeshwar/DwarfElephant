/**
 * This Output is required to couple the DwarfElephant package over a fork
 * interface to the OpenDA framework.
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTSTATEOUTPUT_H
#define DWARFELEPHANTSTATEOUTPUT_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "FileOutput.h"

///-------------------------------------------------------------------------
// Forward declerations
class DwarfElephantStateOutput;

///----------------------------INPUT PARAMETERS-----------------------------
template <>
InputParameters validParams<DwarfElephantStateOutput>();

///-------------------------------------------------------------------------
class DwarfElephantStateOutput : public FileOutput
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantStateOutput(const InputParameters & parameters);

  /*Methods*/
  virtual void output(const ExecFlagType & type) override;

  std::string filename() override;

//--------------------------------PROTECTED---------------------------------
protected:
  /*Attributes*/
  bool _use_rb;

  std::string _system_name;
  VariableName _var_name;
};

///-------------------------------------------------------------------------
#endif /* DWARFELEPHANTSTATEOUTPUT_H */
