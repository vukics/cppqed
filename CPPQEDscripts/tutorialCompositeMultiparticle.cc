// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"

#include "ParticleCavity.h"
#include "IdenticalParticles.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p);
  particle::ParsPumped pp(p);
  mode::ParsLossy pm(p);
  particlecavity::ParsOrthogonal ppc(p);

  update(p,argc,argv,"--");

  const auto mode{std::make_shared<LossyMode<>>(pm)}; // Free0
  const auto part{std::make_shared<PumpedParticle>(pp)}; // Free1,2 – only one instant

  const auto act{std::make_shared<ParticleOrthogonalToCavity>(mode,part,ppc)}; // only one instant

  const auto psi{std::make_shared<quantumdata::StateVector<3>>(init(pm)*hoState(pp)*hoState(pp))};

  evolve(psi,
         composite::make(
                         _<0,1>(act),
                         _<0,2>(act),
                         _<1,2>(std::make_shared<IdenticalParticles<2>>(part /*,...*/))
                         ),
         pe);
  
}
