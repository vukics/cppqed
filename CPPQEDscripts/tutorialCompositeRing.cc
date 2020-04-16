// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"

#include "ParticleTwoModes.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p);
  particle::Pars pp(p);
  mode::ParsLossy       pmP(p,"P");
  mode::ParsPumpedLossy pmM(p,"M");
  particlecavity::ParsAlong ppcP(p,"P");
  particlecavity::ParsAlong ppcM(p,"M");

  ppcP.modeCav=MFT_PLUS; ppcM.modeCav=MFT_MINUS; 

  auto qmp=updateWithPicture(p,argc,argv,"--");

  particle::Ptr part (make(pp ,qmp));
  mode    ::Ptr plus (make(pmP,qmp));
  mode    ::Ptr minus(make(pmM,qmp));

  evolve<0>(wavePacket(pp)*init(pmP)*init(pmM),
            composite::make(
                            _<1,0>  (std::make_shared<ParticleAlongCavity>(plus ,part,ppcP)),
                            _<2,0>  (std::make_shared<ParticleAlongCavity>(minus,part,ppcM)),
                            _<1,2,0>(std::make_shared<ParticleTwoModes>(plus,minus,part,ppcP,ppcM))
                            ),
            pe);

  {
  const composite::result_of::Make<_<1,0>,_<2,0>,_<1,2,0> >::type
    system=composite::make(
                           _<1,0>  (std::make_shared<ParticleAlongCavity>(plus ,part,ppcP)),
                           _<2,0>  (std::make_shared<ParticleAlongCavity>(minus,part,ppcM)),
                           _<1,2,0>(std::make_shared<ParticleTwoModes>(plus,minus,part,ppcP,ppcM)));

  system->displayParameters(std::cout);
  }
  {
  const auto
    system=composite::make(
                           _<1,0>  (std::make_shared<ParticleAlongCavity>(plus ,part,ppcP)),
                           _<2,0>  (std::make_shared<ParticleAlongCavity>(minus,part,ppcM)),
                           _<1,2,0>(std::make_shared<ParticleTwoModes>(plus,minus,part,ppcP,ppcM)));

  system->displayParameters(std::cout);
  }

}
