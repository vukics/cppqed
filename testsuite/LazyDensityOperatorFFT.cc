#include "Particle.h"
#include "TMP_Tools.h"
#include "impl/DensityOperator.tcc"
#include "impl/LazyDensityOperator.tcc"
#include "impl/LazyDensityOperatorFFT.tcc"

#include <iostream>

#define BOOST_TEST_MODULE LazyDensityOperatorFFT test
#include <boost/test/unit_test.hpp>

typedef quantumdata::StateVector<3>         SV3;
typedef quantumdata::DensityOperator<3>     DO3;
typedef quantumdata::LazyDensityOperator<3> LDO3;

const double eps=1e-12;

BOOST_AUTO_TEST_CASE( LAZY_DENSITY_OPERATOR_FFTRANSFORM )
{
  using tmptools::Vector;
  using namespace particle;
  using namespace quantumdata;
  using namespace blitz;
  using mathutils::fcmp;
  double s1=1.;
  double s2=2.;
  double s3=3.;
  Spatial spacial(3);
  SV3 psi1 = wavePacket(InitialCondition(0.5,0,s1,true),spacial,true)*wavePacket(InitialCondition(0,0,s2,true),spacial,true)*wavePacket(InitialCondition(-0.5,0,s3,true),spacial,true);
  SV3 psi2 = wavePacket(InitialCondition(0.5,0,s1,true),spacial,false)*wavePacket(InitialCondition(0,0,s2,true),spacial,true)*wavePacket(InitialCondition(-0.5,0,s3,true),spacial,false);
  LDO3::Ptr psiX = ffTransform<Vector<0,2> >(psi1,fft::DIR_KX);
  const SV3 *psi3 = dynamic_cast<const SV3 *>(psiX.get());
  BOOST_CHECK(max(abs(psi2()-(*psi3)()))<eps);
  
  DO3 rho1 = dyad(psi1,psi1);
  DO3 rho2 = dyad(psi2,psi2);
  LDO3::Ptr rhoX = ffTransform<Vector<0,2> >(rho1,fft::DIR_KX);
  const DO3 *rho3 = dynamic_cast<const DO3 *>(rhoX.get());
  BOOST_CHECK(max(abs(rho2()-(*rho3)()))<eps);
}