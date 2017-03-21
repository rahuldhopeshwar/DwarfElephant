 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantOfflineStage.h"

template<>
InputParameters validParams<DwarfElephantOfflineStage>()
{
  InputParameters params = validParams<NodalUserObject>();

  params.addParam<bool>("use_displaced", false, "Enable/disable the use of the displaced mesh for the data retrieving.");
  params.addRequiredParam<bool>("store_basis_functions","Determines whether the basis functions are stored or not.");
  params.addParam<bool>("compliant", true, "Determines whether F is equal to the output vector or not.");
  params.addParam<bool>("skip_matrix_assembly_in_rb_system", true, "Determines whether the matrix is assembled in the RB System or in the nl0 system.");
  params.addParam<bool>("skip_vector_assembly_in_rb_system", true, "Determines whether the vectors are assembled in the RB System or in the nl0 system.");
  params.addParam<std::string>("system","nl0","The name of the system that should be read in.");
  params.addRequiredParam<UserObjectName>("initial_rb_userobject", "Name of the UserObject for initializing the RB system.");

  return params;
}

DwarfElephantOfflineStage::DwarfElephantOfflineStage(const InputParameters & params):
  NodalUserObject(params),
  _use_displaced(getParam<bool>("use_displaced")),
  _store_basis_functions(getParam<bool>("store_basis_functions")),
  _skip_matrix_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _skip_vector_assembly_in_rb_system(getParam<bool>("skip_matrix_assembly_in_rb_system")),
  _compliant(getParam<bool>("compliant")),
  _system_name(getParam<std::string>("system")),
  _block_ids(this->blockIDs()),
  _es(_use_displaced ? _fe_problem.getDisplacedProblem()->es() : _fe_problem.es()),
  _sys(_es.get_system<TransientNonlinearImplicitSystem>(_system_name)),
  _initialize_rb_system(getUserObject<DwarfElephantInitializeRBSystem>("initial_rb_userobject")),
  _mesh_ptr(&_fe_problem.mesh())
{
}

void
DwarfElephantOfflineStage::transferAffineOperators(bool _skip_matrix_assembly_in_rb_system, bool _skip_vector_assembly_in_rb_system)
{
  // Transfer the vectors
  if (_skip_vector_assembly_in_rb_system)
  {
    // Transfer the data for the F vectors.
    for(unsigned int _q=0; _q<1; _q++)
      _initialize_rb_system._residuals[_q]->operator=(_sys.get_vector("Re_non_time"));

    // Transfer the data for the output vectors.
    if (_compliant)
    {
      for(unsigned int _q=0; _q<1; _q++)
        _initialize_rb_system._outputs[_q]->operator=(_sys.get_vector("Re_non_time"));
    }
    else if (!_compliant)
      mooseError("Currently, the code handles the compliant case, only.");
  }

  if (_skip_matrix_assembly_in_rb_system)
  {
    // The stiffness matrices are transfered in the RBKernel class.

    // Transfer the inner product matrix
    _initialize_rb_system._rb_con_ptr->get_inner_product_matrix()->close();
    _initialize_rb_system._rb_con_ptr->get_inner_product_matrix()->add(1,*_sys.matrix);
  }
}

void
DwarfElephantOfflineStage::offlineStage()
{
  // This method performs the offline stage of the RB problem.

  // Computation of the reduced basis space.
//  _initialize_rb_system._rb_con_ptr->print_info();
  _initialize_rb_system._rb_con_ptr->train_reduced_basis();

  // Write the offline data to file (xdr format).
//  _rb_con_ptr->get_rb_evaluation().legacy_write_offline_data_to_files();
//
//  // If desired, store the basis functions (xdr format).
//  if (_store_basis_functions)
//  {
//    _rb_con_ptr->get_rb_evaluation().write_out_basis_functions(*_rb_con_ptr);
//  }
//
//  _rb_con_ptr->print_basis_function_orthogonality();
}

void
DwarfElephantOfflineStage::initialize()
{
}

void
DwarfElephantOfflineStage::execute()
{
}

void
DwarfElephantOfflineStage::threadJoin(const UserObject & y)
{

}


void
DwarfElephantOfflineStage::finalize()
{
  transferAffineOperators(_skip_matrix_assembly_in_rb_system, _skip_vector_assembly_in_rb_system);
  offlineStage();
}
