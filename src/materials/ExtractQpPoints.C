#include "ExtractQpPoints.h"

//----------------------------INPUT PARAMETERS-----------------------------
template<>
InputParameters validParams<ExtractQpPoints>()
{
  InputParameters params = validParams<Material>();
  return params;
}

//-------------------------------CONSTRUCTOR-------------------------------
ExtractQpPoints::ExtractQpPoints(const InputParameters & parameters) :
    Material(parameters)
{
}

//-------------------------------------------------------------------------
void
ExtractQpPoints::computeQpProperties()
{
  // NonlinearSystemBase & _nl_sys = _fe_problem.getNonlinearSystemBase();

  std::ofstream _qp_file;
  _qp_file.open("QpPoints.txt", std::ios::app);
}
