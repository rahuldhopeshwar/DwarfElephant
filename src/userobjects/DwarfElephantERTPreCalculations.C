 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantERTPreCalculations.h"

template<>
InputParameters validParams<DwarfElephantERTPreCalculations>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<DenseVector<Real>>("position_A_electrode", "Positions of the first current electrode");
  params.addRequiredParam<DenseVector<Real>>("position_B_electrode", "Positions of the second current electrode");
  params.addRequiredParam<DenseVector<Real>>("position_M_electrode", "Positions of the first potential electrode");
  params.addRequiredParam<DenseVector<Real>>("position_N_electrode", "Positions of the second potential electrode");

  return params;
}

DwarfElephantERTPreCalculations::DwarfElephantERTPreCalculations(const InputParameters & params):
  GeneralUserObject(params),
  _exec_flags(this->execFlags()),
  _position_A_electrode(getParam<DenseVector<Real>>("position_A_electrode")),
  _position_B_electrode(getParam<DenseVector<Real>>("position_B_electrode")),
  _position_M_electrode(getParam<DenseVector<Real>>("position_M_electrode")),
  _position_N_electrode(getParam<DenseVector<Real>>("position_N_electrode"))
{
}

void
DwarfElephantERTPreCalculations::initialize()
{
 if (_exec_flags[0]!=EXEC_INITIAL)
  mooseError("UserObject has to be performed on initial stage.");
}

void
DwarfElephantERTPreCalculations::execute()
{
}

void
DwarfElephantERTPreCalculations::finalize()
{
}

DenseVector<Real> &
DwarfElephantERTPreCalculations::computeGemeotricFactor()
{
  (computeDistance(_position_A_electrode, _position_M_electrode)-=computeDistance(_position_A_electrode, _position_N_electrode)-=
   computeDistance(_position_B_electrode, _position_M_electrode)+= computeDistance(_position_B_electrode, _position_N_electrode))*=(1./libMesh::pi);
  return _geometric_factor;
}

DenseVector<Real> &
DwarfElephantERTPreCalculations::computeDistance(DenseVector<Real> & vec1, DenseVector<Real> & vec2)
{
  return vec2-=vec1;
}
