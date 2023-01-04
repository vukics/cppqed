// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionBinary.h"

#include "ParticleCavity.h"


using namespace std;
using namespace cppqedutils;


typedef quantumdata::StateVector    <2> StateVector    ;
typedef quantumdata::DensityOperator<2> DensityOperator;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  int& conf=p.add("1p1mconf","System configuration code for 1particle1mode",1);

  evolution::Pars<> pe(p); // Driver Parameters
  mode::ParsDrivenDissipative pplm(p); // (Driven) Cavity
  particle::ParsDriven ppp(p); // Driven Particle
  particlecavity::ParsAlong ppci(p); // Particle Cavity Interaction

  dcomp& superposition=p.add("superposition","Particle initial condition |0>+|2>",dcomp(0,0));

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  if (conf!=1) pplm.delta-=ppci.uNot/(isComplex(ppci.modeCav) ? 1. : 2.);

  mode::Ptr mode(make(pplm,qmp));


  particle::Ptr particle(make(ppp,qmp));

  particle::PtrDriven drivenparticle(makeDriven(ppp,qmp));

  std::shared_ptr<particlecavity::Base> particlecavityBase;
  
  switch (conf) {
  case 1:
    // Driven particle moving orthogonal to (driven) cavity
    particlecavityBase.reset(new ParticleOrthogonalToCavity(mode,drivenparticle,ppci));
    break;
  case 2:
    // (Driven) particle moving along cavity. 
    // Atomic pump aligned orthogonal to cavity -> doesn't affect atomic motion directly. 
    // Cavity driven + atom scatters from atomic pump into cavity.
    if (!abs(pplm.eta) && !ppp.vClass) {cerr<<"No driving in the system!"<<endl; return 1;}
    if (!ppp.init.getSig()) {cerr<<"No initial spread specified!"<<endl; return 1;}
    particlecavityBase.reset(new ParticleAlongCavity(mode,particle,ppci,ppp.vClass));
    break;
  case 3:
    // Driven particle moving along cavity. 
    // Atomic pump aligned along cavity -> atomic moves in additional classical potential. 
    // Cavity driven + atom scatters from atomic pump into cavity.
    particlecavityBase.reset(new ParticleAlongCavity(mode,drivenparticle,ppci));
    break;
  case 4:
    // Same as case 3, but atom doesn't scatter from atomic pump to cavity.
    // -> Atomic pump merely a classical potential.
    particlecavityBase.reset(new ParticleAlongCavity(mode,drivenparticle,ppci,0));
    break;
  default:
    cerr<<"Configuration not recognized!"<<endl;
    return 1;
  }

  // ****** ****** ****** ****** ****** ******

  quantumdata::StateVector<1> 
    psiPart=
    abs(superposition)>0 
    ? 
    sqrt(1-sqrAbs(superposition))*hoState(0,ppp.init,particle::Spatial(ppp.fin))+superposition*hoState(2,ppp.init,particle::Spatial(ppp.fin)) 
    :
    init(ppp)
    ;

  StateVector psi=init(pplm)*psiPart; psi.renorm();


  evolve(psi,binary::make(particlecavityBase),pe);

}
