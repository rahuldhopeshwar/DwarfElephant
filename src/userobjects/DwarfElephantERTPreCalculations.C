 ///-------------------------------------------------------------------------
// MOOSE includes (DwarfElephant package)
#include "DwarfElephantERTPreCalculations.h"

registerMooseObject("DwarfElephantApp", DwarfElephantERTPreCalculations);

template<>
InputParameters validParams<DwarfElephantERTPreCalculations>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<unsigned int>("n_electrodes", "Number of electrodes");
  params.addRequiredParam<std::vector<Real>>("position_A_electrode", "Positions of the first current electrode");
  params.addRequiredParam<std::vector<Real>>("position_B_electrode", "Positions of the second current electrode");
  params.addRequiredParam<std::vector<Real>>("position_M_electrode", "Positions of the first potential electrode");
  params.addRequiredParam<std::vector<Real>>("position_N_electrode", "Positions of the second potential electrode");

  return params;
}

DwarfElephantERTPreCalculations::DwarfElephantERTPreCalculations(const InputParameters & params):
  GeneralUserObject(params),
  _exec_flags(this->execFlags()),
  _n_electrodes(getParam<unsigned int>("n_electrodes")),
  _position_A_electrode(getParam<std::vector<Real>>("position_A_electrode")),
  _position_B_electrode(getParam<std::vector<Real>>("position_B_electrode")),
  _position_M_electrode(getParam<std::vector<Real>>("position_M_electrode")),
  _position_N_electrode(getParam<std::vector<Real>>("position_N_electrode"))
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
  setUpElectrodeVectors();
  computeGeometricFactor();
}

void
DwarfElephantERTPreCalculations::finalize()
{
}

void
DwarfElephantERTPreCalculations::setUpElectrodeVectors()
{
  _A_electrode = NumericVector<Number>::build(this->comm());
  _A_electrode->init(_n_electrodes-3,_n_electrodes-3);

  _B_electrode = NumericVector<Number>::build(this->comm());
  _B_electrode->init(_n_electrodes-3,_n_electrodes-3);

  _M_electrode = NumericVector<Number>::build(this->comm());
  _M_electrode->init(_n_electrodes-3,_n_electrodes-3);

  _N_electrode = NumericVector<Number>::build(this->comm());
  _N_electrode->init(_n_electrodes-3,_n_electrodes-3);

  _dist = NumericVector<Number>::build(this->comm());
  _dist->init(_n_electrodes-3,_n_electrodes-3);

  _geometric_factor = NumericVector<Number>::build(this->comm());
  _geometric_factor->init(_n_electrodes-3,_n_electrodes-3);
  _geometric_factor->zero();

  for (unsigned int i = 0; i < _n_electrodes-3; i++)
  {
    _A_electrode->set(i, _position_A_electrode[i]);
    _B_electrode->set(i, _position_B_electrode[i]);
    _M_electrode->set(i, _position_M_electrode[i]);
    _N_electrode->set(i, _position_N_electrode[i]);
  }
}

NumericVector<Number> &
DwarfElephantERTPreCalculations::computeGeometricFactor()
{
  *_geometric_factor = computeDistance(*_A_electrode,*_M_electrode);
  *_geometric_factor -= computeDistance(*_A_electrode,*_N_electrode);
  *_geometric_factor -= computeDistance(*_B_electrode,*_M_electrode);
  *_geometric_factor += computeDistance(*_B_electrode,*_N_electrode);
  _geometric_factor->reciprocal();
  *_geometric_factor *= (1./(2*libMesh::pi));

  return *_geometric_factor;
}

NumericVector<Number> &
DwarfElephantERTPreCalculations::computeDistance(NumericVector<Number> & vec1, NumericVector<Number> & vec2)
{
  _dist->zero();
  _dist = vec2.clone();
  *_dist -= vec1;
   _dist->abs();
  return *_dist;
}
