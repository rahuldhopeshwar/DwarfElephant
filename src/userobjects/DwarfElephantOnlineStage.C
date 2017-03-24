 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantOnlineStage.h"

template<>
InputParameters validParams<DwarfElephantOnlineStage>()
{
  InputParameters params = validParams<GeneralUserObject>();

  return params;
}

DwarfElephantOnlineStage::DwarfElephantOnlineStage(const InputParameters & params):
  GeneralUserObject(params)
{
}

void
DwarfElephantOnlineStage::onlineStage()
{
  ////    _rb_eval.legacy_read_offline_data_from_files();
////    RBParameters _online_mu;
////
////    _online_mu.set_value("mu0", _online_mu0_parameters);
////    _rb_eval.set_parameters(_online_mu);
////    _rb_eval.rb_solve(_online_N);
////
////    _rb_eval.print_parameters();
}

void
DwarfElephantOnlineStage::initialize()
{
}

void
DwarfElephantOnlineStage::execute()
{
}

void
DwarfElephantOnlineStage::finalize()
{
}
