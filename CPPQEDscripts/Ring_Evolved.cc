// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution_.h"

#include "Simulated.h"

#include "ParticleTwoModes.h"

using namespace std       ;
using namespace cppqedutils;
using namespace trajectory;

typedef CArray<1> Array;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  evolution::Pars<> pt(p);

  particle::Pars pp(p);
  mode::ParsDrivenDissipative pmP(p,"P");
  mode::ParsDrivenDissipative pmM(p,"M");
  particlecavity::ParsAlong ppcP(p,"P");
  particlecavity::ParsAlong ppcM(p,"M");

  ppcP.modeCav=MFT_PLUS; ppcM.modeCav=MFT_MINUS; 

  update(p,argc,argv,"--");

  Array alpha(2); alpha=pmP.minit,pmM.minit;

  run(simulated::makeBoost(alpha,[&](const Array& b, Array& dbdt, double) {
    const double x=PI*pp.init.getX0();

    const dcomp
      z1=1i*(pmP.delta-ppcP.uNot*sqrAbs(modeFunction(ppcP.modeCav,x)))-pmP.kappa,
      z2=1i*(pmM.delta-ppcM.uNot*sqrAbs(modeFunction(ppcM.modeCav,x)))-pmM.kappa,
      g=sign(ppcP.uNot)*sqrt(ppcP.uNot*ppcM.uNot)*conj(modeFunction(ppcP.modeCav,x))*modeFunction(ppcM.modeCav,x);

    dbdt=
      z1*b(0)+pmP.eta-1i*     g *b(1),
      z2*b(1)+pmM.eta-1i*conj(g)*b(0);
  },{"alphaP","alphaM"},1e-6,pt),pt);

}
