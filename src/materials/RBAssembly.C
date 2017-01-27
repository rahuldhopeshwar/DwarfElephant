#include "RBAssembly.h"
#include "Assembly.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"


template<>
InputParameters validParams<RBAssembly>()
{
  InputParameters params = validParams<Material>();
  return params;
}

RBAssembly::RBAssembly(const InputParameters & parameters) :
    Material(parameters)
{
}

void
RBAssembly::constructJacobianAndResidual()
{
  // Grab reference to linear Lagrange finite element object pointer,
  // currently this is always a linear Lagrange element, so this might need to
  // be generalized if we start working with higher-order elements...
  FEBase * & fe(_assembly.getFE(getParam<bool>("linear_shape_fcns") ? FEType(FIRST, LAGRANGE) : FEType(SECOND, LAGRANGE), _current_elem->dim()));

  // Grab references to FE object's mapping data from the _subproblem's FE object
  const std::vector<Real> & dxidx(fe->get_dxidx());
  const std::vector<Real> & dxidy(fe->get_dxidy());
  const std::vector<Real> & dxidz(fe->get_dxidz());
  const std::vector<Real> & detadx(fe->get_detadx());
  const std::vector<Real> & detady(fe->get_detady());
  const std::vector<Real> & detadz(fe->get_detadz());
  const std::vector<Real> & dzetadx(fe->get_dzetadx());
  const std::vector<Real> & dzetady(fe->get_dzetady());
  const std::vector<Real> & dzetadz(fe->get_dzetadz());

//  const RealTensorValue jacobian = (1/dxidx[qp], 1/dxidy[qp], 1/dxidz[qp],
//                                    1/detadx[qp], 1/detady[qp], 1/detadz[qp],
//                                    1/dzetadx[qp], 1/dzetady[qp], 1/dzetadz[qp]);

  for (unsigned int qp = 0; qp < _qrule->n_points(); qp++)
  {

    // Bounds checking on element data and putting into vector form
    mooseAssert(qp < dxidx.size(), "Insufficient data in dxidx array!");
    mooseAssert(qp < dxidy.size(), "Insufficient data in dxidy array!");
    mooseAssert(qp < dxidz.size(), "Insufficient data in dxidz array!");
    if (_mesh.dimension() >= 2)
    {
      mooseAssert(qp < detadx.size(), "Insufficient data in detadx array!");
      mooseAssert(qp < detady.size(), "Insufficient data in detady array!");
      mooseAssert(qp < detadz.size(), "Insufficient data in detadz array!");
    }
    if (_mesh.dimension() >= 3)
    {
      mooseAssert(qp < dzetadx.size(), "Insufficient data in dzetadx array!");
      mooseAssert(qp < dzetady.size(), "Insufficient data in dzetady array!");
      mooseAssert(qp < dzetadz.size(), "Insufficient data in dzetadz array!");
    }
  }
}

void
RBAssembly::computeProperties()
{
}

