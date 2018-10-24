// MOOSE includes
#include "NonlinearSystemBase.h"
#include "FEProblem.h"

#include "DwarfElephantRBOutput.h"

registerMooseObject("DwarfElephantApp", DwarfElephantRBOutput);

template <>
InputParameters
validParams<DwarfElephantRBOutput>()
{
  InputParameters params = validParams<FileOutput>();

  return params;
}

DwarfElephantRBOutput::DwarfElephantRBOutput(const InputParameters & parameters) :
    FileOutput(parameters),
    _mesh_ptr(&_problem_ptr->mesh())
{
}


void
DwarfElephantRBOutput::output(const ExecFlagType & /*type*/)
{
//   Moose::perf_log.push("write_exodus()", "Execution");
//    DwarfElephantRBEvaluationSteadyState _rb_eval(_mesh_ptr->comm() , *_problem_ptr);
//   _initialize_rb_system._rb_con_ptr->set_rb_evaluation(_rb_eval);

   std::string _systems_for_print[] = {"RBSystem"};
   const std::set<std::string>  _system_names_for_print (_systems_for_print, _systems_for_print+sizeof(_systems_for_print)/sizeof(_systems_for_print[0]));

//   _rb_eval.read_in_basis_functions(*_initialize_rb_system._rb_con_ptr);
//   _initialize_rb_system._rb_con_ptr->load_rb_solution();

   ExodusII_IO(_mesh_ptr->getMesh()).write_equation_systems(_file_base + ".e", *_es_ptr, &_system_names_for_print);

}
