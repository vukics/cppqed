// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"

#include "ParticleCavity.h"
#include "IdenticalParticles.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p);
  particle::ParsDriven pp(p);
  mode::ParsDissipative pm(p);
  particlecavity::ParsOrthogonal ppc(p);

  update(p,argc,argv,"--");

  auto mode{std::make_shared<DissipativeMode<>>(pm)}; // Free0
  auto part{std::make_shared<DrivenParticle>(pp)}; // Free1,2 – only one instant

  auto act{std::make_shared<ParticleOrthogonalToCavity>(mode,part,ppc)}; // only one instant

  evolve(init(pm)*hoState(pp)*hoState(pp),
         composite::make(
                         _<0,1>(act),
                         _<0,2>(act),
                         _<1,2>(std::make_shared<IdenticalParticles<2>>(part /*,...*/))
                         ),
         pe);
  
}
