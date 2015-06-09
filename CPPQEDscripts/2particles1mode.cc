// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"

#include "ParticleCavity.h"
#include "IdenticalParticles.h"


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars pe(p);
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
                         _<0,1>(act),_<0,2>(act),
                         _<1,2>(IdenticalParticles<2>(part))
                         ),
         pe);
}
