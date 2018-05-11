// MOOSE includes
#include "NonlinearSystemBase.h"
#include "FEProblem.h"

#include "DwarfElephantDakotaOutput.h"

template <>
InputParameters
validParams<DwarfElephantDakotaOutput>()
{
  InputParameters params = validParams<FileOutput>();
  params.addParam<std::vector<PostprocessorName>>("postprocessor","Defines the name of the postprocessor(s) you want to use.");
  params.addParam<std::string>("delimiter", "   ", "Defines the delimiter.");
  params.addParam<std::string>("simulation_type", "steady", "Determines whether the simulation is steady state or transient.");
  params.addParam<bool>("use_rb", false, "Defines whether the RB or FE method is used.");
  params.addParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system");
  params.addParam<UserObjectName>("offline_online_rb_userobject", "Name of the UserObject for the  offline and online stage of the RB system");

  return params;
}

DwarfElephantDakotaOutput::DwarfElephantDakotaOutput(const InputParameters & parameters) :
    FileOutput(parameters),
    _delimiter(getParam<std::string>("delimiter")),
    _simulation_type(getParam<std::string>("simulation_type")),
    _use_rb(getParam<bool>("use_rb"))
{
  if(_use_rb)
  {
    _initialize_rb_system_name = getParam<UserObjectName>("initial_rb_userobject");
    _offline_online_rb_system_name = getParam<UserObjectName>("offline_online_rb_userobject");
  } else {
    _postprocessor_name = getParam<std::vector<PostprocessorName>>("postprocessor");
  }
}

// void
// DwarfElephantDakotaOutput::initialSetup(){
//   CSV::initialSetup();
// }

void
DwarfElephantDakotaOutput::output(const ExecFlagType & /*type*/)
{
  Moose::perf_log.push("DakotaOutput()", "Output");
  // This result file enables the use of MOOSE as a forward simulator within Dakota.
  // Which output parameters are printed to the result file can be controlled over the MOOSE input file.

//  if (type == EXEC_TIMESTEP_END)
//  {
  if(processor_id() == 0){
    std::ofstream dakota_file;
    dakota_file.open(filename(), std::ios::app);

    if(!_use_rb)
    {
      for(unsigned int i = 0; i < _postprocessor_name.size(); i++)
        if(i < _postprocessor_name.size()-1)
          dakota_file << _problem_ptr->getPostprocessorValue(_postprocessor_name[i]) << " f"<< _delimiter;
        else
          dakota_file << _problem_ptr->getPostprocessorValue(_postprocessor_name[i]) << " f"<< std::endl;
    } else {
      if(_simulation_type == "steady")
      {
        const DwarfElephantInitializeRBSystemSteadyState & _initialize_rb_system = _problem_ptr->getUserObject<DwarfElephantInitializeRBSystemSteadyState>(_initialize_rb_system_name);
        const DwarfElephantOfflineOnlineStageSteadyState & _offline_online_rb_system = _problem_ptr->getUserObject<DwarfElephantOfflineOnlineStageSteadyState>(_offline_online_rb_system_name);

        for (unsigned int i = 0; i != _initialize_rb_system._n_outputs; i++)
          if(i < _initialize_rb_system._n_outputs-1)
            dakota_file << _offline_online_rb_system._RB_outputs[i] << " f"<< _delimiter;
          else
            dakota_file << _offline_online_rb_system._RB_outputs[i] << " f"<< std::endl;
      } else{
        const DwarfElephantInitializeRBSystemTransient & _initialize_rb_system = _problem_ptr->getUserObject<DwarfElephantInitializeRBSystemTransient>(_initialize_rb_system_name);
        const DwarfElephantOfflineOnlineStageTransient & _offline_online_rb_system = _problem_ptr->getUserObject<DwarfElephantOfflineOnlineStageTransient>(_offline_online_rb_system_name);

        for (unsigned int _time_step = 0; _time_step <= _initialize_rb_system._rb_con_ptr->get_n_time_steps(); _time_step++)
        {
          for (unsigned int i = 0; i < _initialize_rb_system._n_outputs; i++)
          {
            if(i < _initialize_rb_system._n_outputs-1)
              dakota_file << _offline_online_rb_system._RB_outputs_all_timesteps[_time_step][i] << " f"<< _delimiter;
            else
              dakota_file << _offline_online_rb_system._RB_outputs_all_timesteps[_time_step][i] << " f"<< std::endl;
            }
         }
         if (_offline_online_rb_system._output_file)
          mooseError("This Output class using the RB method is specifically designed to work efficiently for repetitive forward simulations. Therefore the output of Exodus files is not desired and the combination of outputting Exodus files and DakotaOutputs is therefore not supported.");
      }
    }
 // }
 // std::string deleteline = "0 f";
 // std::string line;
 //
 // std::ifstream dakota_file_rb;
 // dakota_file_rb.open(_file_path + _result_file_name + ".out", std::ios::app);
 //
 // while (std::getline(dakota_file_rb,line))
 // {
 //   line.replace(line.find(deleteline),deleteline.length(),"");
 //   dakota_file << line << std::endl;
 // }
    dakota_file.close();
  }
  Moose::perf_log.pop("DakotaOutput()", "Output");
}



std::string
DwarfElephantDakotaOutput::filename()
{
  return _file_base;
}
