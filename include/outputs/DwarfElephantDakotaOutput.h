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
// #include "CSV.h"
// #include "DwarfElephantFormattedTable.h"
#include "FileOutput.h"

///-------------------------------------------------------------------------
// Forward declerations
class DwarfElephantDakotaOutput;

///----------------------------INPUT PARAMETERS-----------------------------
template <>
InputParameters validParams<DwarfElephantDakotaOutput>();

///-------------------------------------------------------------------------
// class DwarfElephantDakotaOutput : public CSV
class DwarfElephantDakotaOutput : public FileOutput
{
//----------------------------------PUBLIC----------------------------------
public:
  DwarfElephantDakotaOutput(const InputParameters & parameters);

  /*Methods*/
  virtual void output(const ExecFlagType & type) override;
  // virtual void outputVectorPostprocessors() override;
  // virtual void outputPostprocessors() override;
  virtual std::string filename() override;
  // virtual void initialSetup() override;

//--------------------------------PROTECTED---------------------------------
protected:
  /*Attributes*/
  // bool _write_all_table;
  // bool _write_vector_table;
  // bool _align;
  // unsigned int _precision;
  std::string _delimiter;
  // bool _sort_columns;
  std::vector<PostprocessorName> _postprocessor_name;

};

///-------------------------------------------------------------------------
#endif /* DWARFELEPHANTDAKOTAOUTPUT_H */
