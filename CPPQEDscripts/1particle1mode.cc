#include "EvolutionBinary.h"

#include "ParticleCavity.h"


using namespace std      ;
using namespace mathutils;


typedef quantumdata::StateVector    <2> StateVector    ;
typedef quantumdata::DensityOperator<2> DensityOperator;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  int& conf=p.add("1p1mconf","System configuration code for 1particle1mode",1);

  evolution::Pars pe(p); // Driver Parameters
  mode::ParsPumpedLossy pplm(p); // (Pumped) Cavity
  particle::ParsPumped ppp(p); // Pumped Particle
  particlecavity::ParsAlong ppci(p); // Particle Cavity Interaction

  dcomp& superposition=p.add("superposition","Particle initial condition |0>+|2>",dcomp(0,0));

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  if (conf!=1) pplm.delta-=ppci.uNot/(isComplex(ppci.modeCav) ? 1. : 2.);

  mode::Ptr mode(make(pplm,qmp));


  particle::Ptr particle(make(ppp,qmp));

  particle::PtrPumped pumpedparticle(makePumped(ppp,qmp));

  particlecavity::Base* particlecavityBase;
  switch (conf) {
  case 1:
    // Pumped particle moving orthogonal to (pumped) cavity
    particlecavityBase=new ParticleOrthogonalToCavity(mode,pumpedparticle,ppci);
    break;
  case 2:
    // (Pumped) particle moving along cavity. 
    // Atomic pump aligned orthogonal to cavity -> doesn't affect atomic motion directly. 
    // Cavity pumped + atom scatters from atomic pump into cavity.
    if (!abs(pplm.eta) && !ppp.vClass) {cerr<<"No driving in the system!"<<endl; return 1;}
    if (!ppp.init.getSig()) {cerr<<"No initial spread specified!"<<endl; return 1;}
    particlecavityBase=new ParticleAlongCavity(mode,particle,ppci,ppp.vClass);
    break;
  case 3:
    // Pumped particle moving along cavity. 
    // Atomic pump aligned along cavity -> atomic moves in additional classical potential. 
    // Cavity pumped + atom scatters from atomic pump into cavity.
    particlecavityBase=new ParticleAlongCavity(mode,pumpedparticle,ppci);
    break;
  case 4:
    // Same as case 3, but atom doesn't scatter from atomic pump to cavity.
    // -> Atomic pump merely a classical potential.
    particlecavityBase=new ParticleAlongCavity(mode,pumpedparticle,ppci,0);
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

  StateVector psi=init(pplm)*psiPart;
  
  psi.renorm();


  evolve<0>(psi,binary::make(*particlecavityBase),pe);

}
