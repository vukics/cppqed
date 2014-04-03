// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"

#include "ParticleTwoModes.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars pe(p);
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

  quantumdata::StateVector<3> psi(wavePacket(pp)*init(pmP)*init(pmM));

  evolve<0>(psi,
            composite::make(
                            _<1,0>  (ParticleAlongCavity(plus ,part,ppcP)),
                            _<2,0>  (ParticleAlongCavity(minus,part,ppcM)),
                            _<1,2,0>(ParticleTwoModes(plus,minus,part,ppcP,ppcM))
                            ),
            pe);

  {
  ParticleAlongCavity pac1(plus ,part,ppcP), pac2(minus,part,ppcM);
  ParticleTwoModes ptm(plus,minus,part,ppcP,ppcM);

  const composite::result_of::Make<_<1,0>,_<2,0>,_<1,2,0> >::type
    system=composite::make(
                           _<1,0>  (pac1),
                           _<2,0>  (pac2),
                           _<1,2,0>(ptm));

  system->displayParameters(std::cout);
}
  {
  ParticleAlongCavity pac1(plus ,part,ppcP), pac2(minus,part,ppcM);
  ParticleTwoModes ptm(plus,minus,part,ppcP,ppcM);

  const auto
    system=composite::make(
                           _<1,0>  (pac1),
                           _<2,0>  (pac2),
                           _<1,2,0>(ptm));

  system->displayParameters(std::cout);
}

}
