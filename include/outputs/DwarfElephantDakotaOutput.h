/**
 * This Output is required to couple the DwarfElephant package over a fork
 * interface to the multitool kit DAKOTA (developed by Sandia National
 * Labratories).
 */

///-------------------------------------------------------------------------
#ifndef DWARFELEPHANTDAKOTAOUTPUT_H
#define DWARFELEPHANTDAKOTAOUTPUT_H

///---------------------------------INCLUDES--------------------------------
// MOOSE includes
#include "FileOutput.h"

///-------------------------------------------------------------------------
// Forward declerations
class DwarfElephantDakotaOutput;

///----------------------------INPUT PARAMETERS-----------------------------
template <>
InputParameters validParams<DwarfElephantDakotaOutput>();

///-------------------------------------------------------------------------
class DwarfElephantDakotaOutput : public FileOutput
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantDakotaOutput(const InputParameters & parameters);

  /*Methods*/
  virtual void output(const ExecFlagType & type) override;

//--------------------------------PROTECTED---------------------------------
protected:
  /*Attributes*/
  std::string _result_file_name;
  std::string _file_path;

  PostprocessorName _postprocessor_name;
};

///-------------------------------------------------------------------------
#endif /* DWARFELEPHANTDAKOTAOUTPUT_H */
