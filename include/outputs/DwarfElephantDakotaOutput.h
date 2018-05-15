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
#include "DwarfElephantInitializeRBSystemSteadyState.h"
#include "DwarfElephantInitializeRBSystemTransient.h"
#include "DwarfElephantOfflineOnlineStageSteadyState.h"
#include "DwarfElephantOfflineOnlineStageTransient.h"

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
  virtual std::string filename() override;
  // virtual void initialSetup() override;

//--------------------------------PROTECTED---------------------------------
protected:
  /*Attributes*/
  std::string _delimiter;
  std::string _simulation_type;
  std::vector<PostprocessorName> _postprocessor_name;
  UserObjectName _offline_online_rb_system_name;

  bool _use_rb;

};

///-------------------------------------------------------------------------
#endif /* DWARFELEPHANTDAKOTAOUTPUT_H */
