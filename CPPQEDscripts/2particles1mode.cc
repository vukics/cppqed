#include "EvolutionComposite.h"

#include "ParticleCavity.h"
#include "IdenticalParticles.h"


int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsEvolution pe(p);
  particle::ParsPumped  pp(p);
  mode::ParsPumpedLossy pm(p);
  particlecavity::ParsOrthogonal ppc(p);
  
  update(p,argc,argv,"--");

  PumpedLossyMode<> mode(pm); // Free0
  PumpedParticle    part(pp); // Free1,2 - only one instant

  ParticleOrthogonalToCavity act(mode,part,ppc); // only one instant

  quantumdata::StateVector<3> psi(init(pm)*init(pp)*init(pp));
  
  evolve(psi,
         composite::make(
                         Act<0,1>(act),Act<0,2>(act),
                         Act<1,2>(IdenticalParticles<2>(part))
                         ),
         pe);
}
