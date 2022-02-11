// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"

#include "ParticleCavity.h"
#include "IdenticalParticles.h"


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p);
  particle::ParsPumped  pp(p);
  mode::ParsPumpedLossy pm(p);
  particlecavity::ParsOrthogonal ppc(p);
  
  update(p,argc,argv,"--");

  auto mode{std::make_shared<PumpedLossyMode<>>(pm)}; // Free0
  auto part{std::make_shared<PumpedParticle>(pp)}; // Free1,2 - only one instant

  auto act{std::make_shared<const ParticleOrthogonalToCavity>(mode,part,ppc)}; // only one instant

  evolve(init(pm)*init(pp)*init(pp),
         composite::make(
                         _<0,1>(act),_<0,2>(act),
                         _<1,2>(std::make_shared<IdenticalParticles<2>>(part))
                         ),
         pe);
}
