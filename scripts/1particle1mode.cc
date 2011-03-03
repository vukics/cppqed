#include "EvolutionHigh.h"

#include "Mode.h"
#include "ParsMode.h"
#include "ParsParticle.h"
#include "Particle.h"
#include "ParticleInitialCondition.h"

#include "ParsParticleCavity.h"
#include "ParticleCavity.h"

#include "BinarySystem.h"

#include "ModeFunctionType.h"


using namespace std      ;
using namespace mathutils;


typedef quantumdata::StateVector    <2> StateVector    ;
typedef quantumdata::DensityOperator<2> DensityOperator;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  int& conf=p.add("1p1mconf","System configuration code for 1particle1mode",1);

  ParsEvolution pe(p); // Driver Parameters
  mode::ParsPumpedLossy pplm(p); // (Pumped) Cavity
  particle::ParsPumped ppp(p); // Pumped Particle
  particlecavity::ParsAlong ppci(p); // Particle Cavity Interaction

  QM_Picture& qmp=p.add("picture","QM_Picture for particle (IP=UIP or Sch) and field (IP, UIP, or Sch)",QMP_IP);
  // Should be used only out of curiosity (or testing, of course).

  dcomp& superposition=p.add("superposition","Particle initial condition |0>+|2>",dcomp(0,0));

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER || pe.evol==EM_MASTER_FAST /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  if (conf!=1) pplm.delta-=ppci.uNot/(isComplex(ppci.modeCav) ? 1. : 2.);

  mode::SmartPtr mode(mode::maker(pplm,qmp));

  Particle    particle   (ppp);
  ParticleSch particleSch(ppp);

  PumpedParticle    pumpedparticle   (ppp);
  PumpedParticleSch pumpedparticleSch(ppp);


  ParticleBase& particleBase=(
			      qmp==QMP_SCH 
			      ? 
			      static_cast<ParticleBase&>(particleSch) 
			      : 
			      static_cast<ParticleBase&>(particle   )
			      );

  PumpedParticleBase& pumpedparticleBase=(
					  qmp==QMP_SCH 
					  ? 
					  static_cast<PumpedParticleBase&>(pumpedparticleSch) 
					  : 
					  static_cast<PumpedParticleBase&>(pumpedparticle   )
					  );

  particlecavity::Base* particlecavityBase;
  switch (conf) {
  case 1:
    // Pumped particle moving orthogonal to (pumped) cavity
    particlecavityBase=new ParticleOrthogonalToCavity(*mode.get(),pumpedparticleBase,ppci);
    break;
  case 2:
    // (Pumped) particle moving along cavity. 
    // Atomic pump aligned orthogonal to cavity -> doesn't affect atomic motion directly. 
    // Cavity pumped + atom scatters from atomic pump into cavity.
    if (!abs(pplm.eta) && !ppp.vClass) {cerr<<"No driving in the system!"<<endl; return 1;}
    if (!ppp.init.getSig()) {cerr<<"No initial spread specified!"<<endl; return 1;}
    particlecavityBase=new ParticleAlongCavity(*mode.get(),particleBase,ppci,ppp.vClass);
    break;
  case 3:
    // Pumped particle moving along cavity. 
    // Atomic pump aligned along cavity -> atomic moves in additional classical potential. 
    // Cavity pumped + atom scatters from atomic pump into cavity.
    particlecavityBase=new ParticleAlongCavity(*mode.get(),pumpedparticleBase,ppci);
    break;
  case 4:
    // Same as case 3, but atom doesn't scatter from atomic pump to cavity.
    // -> Atomic pump merely a classical potential.
    particlecavityBase=new ParticleAlongCavity(*mode.get(),pumpedparticleBase,ppci,0);
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


  evolve(psi,BinarySystem(*particlecavityBase),pe,tmptools::Vector<0>());


  /*
  BinarySystem system(*particlecavityBase);

  DensityOperator rho(system.getDimensions(),false); rho()=0;

  for (int i=0; i<1; i++) psi.addTo(rho); rho/=1.;

  // cerr<<negPT(rho,tmptools::Vector<0>())<<endl;

    {
    using namespace structure;

    AveragedCommon::Averages averages(Averaged<2>::average(rho,&system));

    AveragedCommon::process(averages,&system);

    cout<<averages<<endl;
  }

  quantumtrajectory::EnsembleMCWF<2,tmptools::Vector<0> > traj(psi,system,pe);

  traj.display();
  
  */

  } catch (ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}
}