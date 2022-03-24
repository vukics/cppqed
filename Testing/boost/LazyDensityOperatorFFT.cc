// Copyright Raimar Sandner 2012–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Particle.h"
#include "TMP_Tools.h"
#include "DensityOperator.h"
#include "LazyDensityOperator.h"
#include "LazyDensityOperatorFFT.h"

#include <iostream>
#include <boost/test/unit_test.hpp>

typedef quantumdata::StateVector<3>         SV3;
typedef quantumdata::DensityOperator<3>     DO3;
typedef quantumdata::LazyDensityOperator<3> LDO3;


BOOST_AUTO_TEST_CASE( LAZY_DENSITY_OPERATOR_FFTRANSFORM , * boost::unit_test::tolerance(1e-12) )
{
  using tmptools::Vector;
  using namespace particle;
  using namespace quantumdata;
  using namespace blitz;
  using cppqedutils::fcmp;
  InitialCondition ic1(0.5,0,1,true);
  InitialCondition ic2(0.5,0,2,true);
  InitialCondition ic3(0.5,0,3,true);
  Spatial spacial(3);
  SV3 psi1 = wavePacket(ic1,spacial,true)*wavePacket(ic2,spacial,true)*wavePacket(ic3,spacial,true);
  SV3 psi2 = wavePacket(ic1,spacial,false)*wavePacket(ic2,spacial,true)*wavePacket(ic3,spacial,false);
  auto psiX = ffTransform<Vector<0,2> >(psi1,fft::DIR_KX);
  const SV3 *psi3 = dynamic_cast<const SV3 *>(psiX.get());
  BOOST_TEST(max(abs(psi2.getArray()-psi3->getArray())) == 0);
  
  DO3 rho1 = dyad(psi1,psi1);
  DO3 rho2 = dyad(psi2,psi2);
  auto rhoX = ffTransform<Vector<0,2> >(rho1,fft::DIR_KX);
  const DO3 *rho3 = dynamic_cast<const DO3 *>(rhoX.get());
  BOOST_TEST(max(abs(rho2.getArray()-rho3->getArray())) == 0);
}
