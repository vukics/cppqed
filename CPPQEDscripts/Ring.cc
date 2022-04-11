// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"

#include "ParticleTwoModes.h"


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p);
  particle::Pars pp(p);
  mode::ParsPumpedLossy pmP(p,"P");
  mode::ParsPumpedLossy pmM(p,"M");
  particlecavity::ParsAlong ppcP(p,"P");
  particlecavity::ParsAlong ppcM(p,"M");

  ppcP.modeCav=MFT_PLUS; ppcM.modeCav=MFT_MINUS; 

  QM_Picture& qmp=updateWithPicture(p,argc,argv);

  pmP.delta-=ppcP.uNot/(isComplex(ppcP.modeCav) ? 1. : 2.);
  pmM.delta-=ppcM.uNot/(isComplex(ppcM.modeCav) ? 1. : 2.);

  particle::Ptr part (make(pp ,qmp));
  mode    ::Ptr plus (make(pmP,qmp));
  mode    ::Ptr minus(make(pmM,qmp));

  evolve<0>(init(pp)*init(pmP)*init(pmM),
            composite::make(
                            _<1,0>  (std::make_shared<ParticleAlongCavity>(plus ,part,ppcP)),
                            _<2,0>  (std::make_shared<ParticleAlongCavity>(minus,part,ppcM)),
                            _<1,2,0>(std::make_shared<ParticleTwoModes>(plus,minus,part,ppcP,ppcM))
                            ),
            pe);

}
