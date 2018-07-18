#include "ExampleGeneralUserObject.h"

template<>
InputParameters validParams<ExampleGeneralUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addRequiredParam<NonlinearVariableName>("variable","The name of the variable that this UserObject operates on");
  return params;
}

ExampleGeneralUserObject::ExampleGeneralUserObject(const InputParameters & params):
  GeneralUserObject(params),
  MooseVariableInterface<Real>(this, false, "variable"),
  _var(*mooseVariable()),
  _JxW(_assembly.JxW()),
  _coord(_assembly.coordTransformation()),
  _phi(_assembly.phi(_var)),
  _test(_var.phi()),
  _qrule(_assembly.qRule())
{
  std::cout << "ExampleGeneralUserObject created" << std::endl;
}

void ExampleGeneralUserObject::initialize()
{
  std::cout << "Started EIM InnerProductMatrix Calculation" << std::endl;
  
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();
  //_inner_product_matrix->init(ke.m(), ke.n(), ke.m(), ke.n());

  for (_i = 0; _i < _phi.size(); _i++)
    for (_j = 0; _j < _phi.size(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _local_ke(_i,_j) += _JxW[_qp] * _coord[_qp] * _phi[_i][_qp] * _phi[_j][_qp];

  //_inner_product_matrix -> add_matrix(_local_ke,_var.dofIndices());
  std::cout << "EIM InnerProductMatrix calculation complete" << std::endl;
  std::cout << "_local_ke matrix calculated using _phi" << std::endl;
  _local_ke.print();


  _local_ke.zero();
  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _test.size(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _local_ke(_i,_j) += _JxW[_qp] * _coord[_qp] * _test[_i][_qp] * _test[_j][_qp];

  std::cout << "_local_ke matrix calculated using _test" << "_test[" << 0 << "][" << 0 << "] = " << _test[0][0] << std::endl;
  _local_ke.print();
}

void ExampleGeneralUserObject::execute()
{

}

void ExampleGeneralUserObject::finalize()
{

}

