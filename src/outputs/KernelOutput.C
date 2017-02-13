///---------------------------------INCLUDE---------------------------------
// Moose includes
#include "FEProblem.h"
#include "FileMesh.h"
#include "MooseApp.h"

// Moose includes (DwarfElephant package)
#include "KernelOutput.h"

///----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<KernelOutput>()
{

  // Get the parameters from the parent object
  InputParameters params = validParams<AdvancedOutput< FileOutput > >();
  params += AdvancedOutput< FileOutput >::enableOutputTypes("scalar postprocessor input system_information");

  params.addRequiredParam<AuxVariableName>("variable", "The name of \
                                                  the variable that this \
                                                  class should output.");

  //params.set<MultiMooseEnum>("execute_on", /*quiet_mode=*/true).push_back("initial timestep_begin linear nonlinear failed");
  params.set<MultiMooseEnum>("execute_system_information_on") = "initial timestep_end";
  params.set<MultiMooseEnum>("execute_postprocessors_on") = "initial timestep_end";
  params.set<MultiMooseEnum>("execute_scalars_on") = "initial timestep_end";

  params.set<std::string>("built_by_action") = "add_user_object";


  return params;
}

///-------------------------------------------------------------------------
KernelOutput::KernelOutput(const InputParameters & parameters) :
    AdvancedOutput<FileOutput>(parameters),
    _tid(parameters.get<THREAD_ID>("_tid")),
//    _non_sys_ptr(&_problem_ptr->getNonlinearSystem()),
    _aux_sys_ptr(&_problem_ptr->getAuxiliarySystem()),
    _n_aux_var(_aux_sys_ptr->nVariables()),
    _var(_problem_ptr->getVariable(_tid, parameters.get<AuxVariableName>("variable"))),
    _nodal(_var.isNodal()),
    _u(_nodal ? _var.nodalSln() : _var.sln()),
    _residual(_var.sys().solution())
//    _var(_aux_sys_ptr->getVariable(_tid, 1))

{
}

///-------------------------------------------------------------------------
void
KernelOutput::output(const ExecFlagType & type)
{
  if (type==EXEC_TIMESTEP_END)
  {
    std::ofstream _file_output;
    _file_output.open("loadVectorF0.xdr");

//    const NumericVector< Number > *& _residual_current = _var.sys().currentSolution();
//    NumericVector< Number >  & _residual = _var.sys().solution();
//    NumericVector< Number > & _residual_old = _var.sys().solutionOld();

//    _file_output << _residual_current << std::endl;
//    _file_output << std::endl;
    _file_output << _residual << std::endl;
//    _file_output << _residual_old << std::endl;

//    _file_output << std::endl;
    _file_output.close();

    _console << std::endl;
    _console << _var.name() << std::endl;
//    _console << _var.sln() << std::endl;
    _console << std::endl;
//    _console << _residual << std::endl;
//    _console << std::endl;

  }
}
